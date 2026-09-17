"""Step 2 — download one light curve per target star and write a manifest.

Usage:
    uv run python -m src.fetch                 # sector from .env (SECTOR=104)
    uv run python -m src.fetch --sector 104 --limit 20   # quick test on 20 stars
    uv run python -m src.fetch --force         # re-download stars already in the manifest

Input : data/targets/sector_{N}.parquet          (from src.targets)
Output: data/lc_raw/**/*.fits                    (lightkurve's mastDownload tree)
        data/targets/sector_{N}_manifest.parquet (tic_id, author, cadence_s, n_points, file_path, status)
"""

from __future__ import annotations

import argparse
import logging
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import dotenv
import joblib
import lightkurve as lk
import pandas as pd
import numpy as np

log = logging.getLogger("fetch")

dotenv.load_dotenv()
DATA_DIR = Path(
    os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data")
)
RAW_LC_DIR = DATA_DIR / "lc_raw"

# Preference order (author, exptime_s). SPOC 2-min beats SPOC 20-s for a transit search
# (same systematics correction, 6x fewer points); TESS-SPOC (PDCSAP) beats QLP (KSPSAP).
PREFERENCE = [
    ("SPOC", 120),
    ("SPOC", 20),
    ("TESS-SPOC", 600),
    ("TESS-SPOC", 200),
    ("QLP", 600),
    ("QLP", 200),
]
AUTHORS = ["SPOC", "TESS-SPOC", "QLP"]
N_RETRIES = 3


@dataclass
class ManifestRow:
    tic_id: int
    author: str = ""
    cadence_s: float = float("nan")
    n_points: int = 0
    file_path: str = ""
    status: str = "error"  # ok | not_found | error
    error: str = ""


def pick_best(search_result: lk.SearchResult) -> lk.SearchResult | None:
    """Return the single best row of a SearchResult according to PREFERENCE, or None."""
    tab = search_result.table
    for author, exptime in PREFERENCE:
        mask: np.ndarray = (tab["author"] == author) & (tab["exptime"] == exptime)
        if mask.any():
            return search_result[int(mask.argmax())]
    return (
        search_result[0] if len(search_result) else None
    )  # unknown cadence: take the first


def fetch_one(tic_id: int, sector: int) -> ManifestRow:
    """Search MAST for `tic_id` in `sector`, download the best product, describe it."""
    row = ManifestRow(tic_id=tic_id)
    for attempt in range(1, N_RETRIES + 1):
        try:
            sr = lk.search_lightcurve(f"TIC {tic_id}", sector=sector, author=AUTHORS)
            best = pick_best(sr)
            if best is None:
                row.status = "not_found"
                return row
            # quality_bitmask="default" drops cadences flagged by the pipeline (momentum dumps etc.)
            lc = best.download(download_dir=str(RAW_LC_DIR), quality_bitmask="default")
            if lc is None:  # lightkurve returns None on a corrupt/unavailable file
                raise RuntimeError("download returned None")
            row.author = str(lc.meta.get("AUTHOR", best.author[0]))
            row.cadence_s = float(best.exptime[0].value)
            row.n_points = int(len(lc.remove_nans()))
            row.file_path = str(lc.meta["FILENAME"])
            row.status = "ok"
            row.error = ""
            return row
        except Exception as e:  # network / MAST hiccups are common; retry with back-off
            row.error = f"{type(e).__name__}: {e}"
            log.warning(
                "TIC %d attempt %d/%d failed: %s", tic_id, attempt, N_RETRIES, row.error
            )
            time.sleep(5 * attempt)
    row.status = "error"
    return row


def fetch_sector(
    sector: int, n_jobs: int = 4, limit: int | None = None, force: bool = False
) -> pd.DataFrame:
    targets = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}.parquet")
    manifest_path = DATA_DIR / "targets" / f"sector_{sector}_manifest.parquet"
    tic_ids = targets["tic_id"].tolist()[:limit]

    # Resume: keep rows already fetched (ok / not_found), redo errors and new stars.
    old = (
        pd.read_parquet(manifest_path)
        if manifest_path.exists() and not force
        else pd.DataFrame(columns=["tic_id", "status"])
    )
    done = set(old.loc[old["status"].isin(["ok", "not_found"]), "tic_id"])
    todo = [t for t in tic_ids if t not in done]
    log.info(
        "sector %d: %d targets, %d already done, %d to fetch",
        sector,
        len(tic_ids),
        len(done),
        len(todo),
    )

    RAW_LC_DIR.mkdir(parents=True, exist_ok=True)
    # n_jobs <= 4: MAST rate-limits; more workers just produce errors.
    rows = joblib.Parallel(n_jobs=n_jobs, verbose=5)(
        joblib.delayed(fetch_one)(t, sector) for t in todo
    )
    new = pd.DataFrame([asdict(r) for r in rows])

    manifest = (
        pd.concat([old[~old["tic_id"].isin(new["tic_id"])], new], ignore_index=True)
        if len(new)
        else old
    )
    manifest = manifest.sort_values("tic_id").reset_index(drop=True)
    manifest.to_parquet(manifest_path, index=False)
    log.info(
        "wrote %s\n%s", manifest_path, manifest["status"].value_counts().to_string()
    )
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--sector", type=int, default=int(os.environ.get("SECTOR", 0)) or None
    )
    p.add_argument("--n-jobs", type=int, default=4)
    p.add_argument(
        "--limit", type=int, default=None, help="only the first N targets (for testing)"
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="ignore the existing manifest and redo everything",
    )
    a = p.parse_args()
    if a.sector is None:
        p.error("no --sector given and SECTOR not set in .env")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    logging.getLogger("lightkurve").setLevel(logging.ERROR)  # silence per-file warnings
    m = fetch_sector(a.sector, n_jobs=a.n_jobs, limit=a.limit, force=a.force)
    print(m.head(10).to_string())


if __name__ == "__main__":
    main()
