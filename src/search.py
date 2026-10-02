"""Step 4 — periodic transit search with Transit Least Squares (TLS) on the flattened light curves.

Usage:
    uv run python -m src.search                        # sector from .env
    uv run python -m src.search --sector 104 --limit 50 --n-jobs 8
    uv run python -m src.search --period-min 0.5 --oversampling 5
    uv run python -m src.search --force
    uv run python -m src.search --second-pass          # after a first pass exists
    uv run python -m src.search --second-pass --only-candidates   # only stars whose pass-1 signal passed vetting

Input : data/targets/sector_{N}_detrend.parquet (rows with status == "ok", from src.detrend)
        data/targets/sector_{N}.parquet           (rad, mass per star, from src.targets)
Output: data/search/sector_{N}.parquet   one row per star with the best TLS signal.
        NO threshold is applied here — every star gets a row; thresholds live in src.vet.

Second pass (--second-pass): reads the first-pass table, masks every predicted transit of
each star's best signal (transit_mask, 2x duration) and runs TLS again on the rest, to find a
second planet. Writes a table with the SAME columns to data/search/pass2/sector_{N}_pass2.parquet.

Walkthrough of every parameter and output field: notebooks/pipeline-4-search-walkthrough.ipynb
"""

from __future__ import annotations

import argparse
import logging
import os
import time
import warnings
from dataclasses import asdict, dataclass, field
from pathlib import Path

import dotenv
import joblib
import numpy as np
import pandas as pd
from transitleastsquares import cleaned_array, transit_mask, transitleastsquares

log = logging.getLogger("search")

dotenv.load_dotenv()
DATA_DIR = Path(
    os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data")
)
SEARCH_DIR = DATA_DIR / "search"

# Shortest trial period. The plan says 0.5 d, but ultra-short-period planets around M dwarfs are a
# real population (LHS 3844 b at 0.463 d is in sector 104) and short periods add few grid points.
DEFAULT_PERIOD_MIN = 0.2


@dataclass
class SearchRow:
    tic_id: int
    n_points: int = 0
    time_span: float = np.nan  # days, per star (gaps and edge trimming differ)
    # --- TLS best signal (pass 1) ---
    SDE: float = np.nan
    SDE_raw: float = np.nan
    FAP: float = np.nan
    snr: float = np.nan
    period: float = np.nan  # days
    period_unc: float = np.nan  # inf when the peak sits at the grid edge
    period_min: float = np.nan  # grid actually used (for the vet.py edge test)
    period_max: float = np.nan
    T0: float = np.nan  # BTJD
    duration_h: float = np.nan
    depth_ppm: float = np.nan  # (1 - depth) * 1e6; TLS's raw depth (~0.99) is easy to misread
    depth_mean_ppm: float = np.nan
    rp_rs: float = np.nan
    transit_count: int = 0  # predicted transits in the baseline
    distinct_transit_count: int = 0  # ... of which have data. vet.py must require >= 2
    empty_transit_count: int = 0
    odd_even_mismatch: float = np.nan  # sigma
    n_periods: int = 0
    transit_times: list[float] = field(default_factory=list)
    transit_depths: list[float] = field(default_factory=list)
    # --- bookkeeping ---
    runtime_s: float = np.nan
    status: str = "error"  # ok | error
    error: str = ""


def _worker_logging() -> None:
    """joblib (loky) workers are fresh processes with no logging config — set it up per worker."""
    if not logging.getLogger().handlers:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s [%(processName)s] %(message)s",
        )


def tls_params(
    rad: float,
    mass: float,
    time_span: float,
    period_min: float = DEFAULT_PERIOD_MIN,
    period_max: float | None = None,
    oversampling_factor: int = 3,
) -> dict:
    """power() keyword arguments for one star. rad/mass in solar units, time_span in days."""
    return dict(
        R_star=float(rad),
        M_star=float(mass),
        # bounds must bracket R_star/M_star or TLS raises ValueError; TIC dwarfs go down to 0.128 R, 0.10 M
        R_star_min=0.1,
        R_star_max=1.0,
        M_star_min=0.05,
        M_star_max=1.0,
        period_min=period_min,
        # at least two transits must fit in the baseline; a single dip is not periodic evidence
        period_max=time_span / 2 if period_max is None else period_max,
        oversampling_factor=oversampling_factor,
        duration_grid_step=1.1,
        use_threads=1,  # one thread per star; parallelism goes across stars
        show_progress_bar=False,
        verbose=False,
    )


def _summarize(r) -> dict:
    """Scalars from a TLS result, in the units we store (depth in ppm, duration in hours)."""
    return {
        "SDE": float(r.SDE),
        "SDE_raw": float(r.SDE_raw),
        "FAP": float(r.FAP),
        "snr": float(r.snr),
        "period": float(r.period),
        "period_unc": float(r.period_uncertainty),
        "T0": float(r.T0),
        "duration_h": float(r.duration) * 24,
        "depth_ppm": (1 - float(r.depth)) * 1e6,
        "depth_mean_ppm": (1 - float(r.depth_mean[0])) * 1e6,
        "rp_rs": float(r.rp_rs),
        "transit_count": int(r.transit_count),
        "distinct_transit_count": int(r.distinct_transit_count),
        "empty_transit_count": int(r.empty_transit_count),
        "odd_even_mismatch": float(r.odd_even_mismatch),
        "n_periods": int(len(r.periods)),
    }


def search_one(
    tic_id: int,
    flat_file: str,
    rad: float,
    mass: float,
    period_min: float = DEFAULT_PERIOD_MIN,
    period_max: float | None = None,
    oversampling_factor: int = 3,
    mask: tuple[float, float, float] | None = None,
) -> SearchRow:
    """Run TLS on one flattened light curve and return one row of scalars. Never raises.

    mask = (period [d], duration [d], T0 [BTJD]) of a signal to remove first (second pass).
    """
    _worker_logging()
    row = SearchRow(tic_id=tic_id)
    t0 = time.time()
    try:
        lc = pd.read_parquet(flat_file, columns=["time", "flux", "flux_err"])
        # safety net: drops non-finite entries, casts to float64
        t, y, dy = cleaned_array(lc["time"].values, lc["flux"].values, lc["flux_err"].values)
        if mask is not None:
            # drop every predicted transit of the first-pass signal (2x duration for margin)
            in_transit = transit_mask(t, mask[0], 2 * mask[1], mask[2])
            t, y, dy = t[~in_transit], y[~in_transit], dy[~in_transit]
        row.n_points = len(t)
        row.time_span = float(t.max() - t.min())
        params = tls_params(rad, mass, row.time_span, period_min, period_max, oversampling_factor)
        row.period_min = params["period_min"]
        row.period_max = params["period_max"]

        with warnings.catch_warnings():
            # "N of M transits without data" — we store empty_transit_count instead of printing it
            warnings.simplefilter("ignore", UserWarning)
            r = transitleastsquares(t, y, dy).power(**params)
            for k, v in _summarize(r).items():
                setattr(row, k, v)
            row.transit_times = [float(x) for x in r.transit_times]
            row.transit_depths = [float(x) for x in r.transit_depths]

        row.runtime_s = time.time() - t0
        row.status = "ok"
        log.info(
            "TIC %d: SDE %.1f, P %.4f d, depth %.0f ppm, %d/%d transits with data, %.0f s",
            tic_id,
            row.SDE,
            row.period,
            row.depth_ppm,
            row.distinct_transit_count,
            row.transit_count,
            row.runtime_s,
        )
    except Exception as e:
        row.runtime_s = time.time() - t0
        row.error = f"{type(e).__name__}: {e}"
        log.error("TIC %d failed: %s", tic_id, row.error)
    return row


def output_path(sector: int, second_pass: bool = False) -> Path:
    if second_pass:
        return SEARCH_DIR / "pass2" / f"sector_{sector}_pass2.parquet"
    return SEARCH_DIR / f"sector_{sector}.parquet"


def search_sector(
    sector: int,
    n_jobs: int = 8,
    limit: int | None = None,
    force: bool = False,
    period_min: float = DEFAULT_PERIOD_MIN,
    period_max: float | None = None,
    oversampling_factor: int = 3,
    second_pass: bool = False,
    only_candidates: bool = False,
) -> pd.DataFrame:
    detrend = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}_detrend.parquet")
    targets = pd.read_parquet(
        DATA_DIR / "targets" / f"sector_{sector}.parquet", columns=["tic_id", "rad", "mass"]
    ).set_index("tic_id")
    ok = detrend[detrend["status"] == "ok"][["tic_id", "file_path"]]
    if second_pass:
        # the stars to search are those with a first-pass signal; carry its ephemeris to mask it
        pass1 = pd.read_parquet(output_path(sector))
        pass1 = pass1[pass1["status"] == "ok"][["tic_id", "period", "duration_h", "T0"]]
        ok = ok.merge(pass1, on="tic_id")
        if only_candidates:
            # a second planet is most interesting (and the search cheapest) where the first signal survived vetting
            vetted = pd.read_parquet(SEARCH_DIR / f"sector_{sector}_vetted.parquet", columns=["tic_id", "passes_all"])
            ok = ok[ok["tic_id"].isin(vetted.loc[vetted["passes_all"], "tic_id"])]
    ok = ok.head(limit).join(targets, on="tic_id")
    missing = ok["rad"].isna() | ok["mass"].isna()
    if missing.any():
        log.warning("%d stars have no rad/mass in the target table; using 0.3/0.3", missing.sum())
        ok.loc[missing, ["rad", "mass"]] = ok.loc[missing, ["rad", "mass"]].fillna(0.3)

    out_path = output_path(sector, second_pass)
    old = pd.read_parquet(out_path) if out_path.exists() else None
    if old is not None and not force:
        ok = ok[~ok["tic_id"].isin(old["tic_id"])]
    log.info(
        "sector %d%s: %d flattened light curves, %d to search (period %.2f-%s d, oversampling %d, "
        "%d workers)",
        sector,
        " second pass" if second_pass else "",
        (detrend["status"] == "ok").sum(),
        len(ok),
        period_min,
        "span/2" if period_max is None else f"{period_max:.2f}",
        oversampling_factor,
        n_jobs,
    )

    rows = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(search_one)(
            int(r.tic_id), r.file_path, float(r.rad), float(r.mass),
            period_min, period_max, oversampling_factor,
            (float(r.period), float(r.duration_h) / 24, float(r.T0)) if second_pass else None,
        )
        for r in ok.itertuples()
    )
    summary = pd.DataFrame([asdict(r) for r in rows])
    if old is not None and not force:
        summary = pd.concat(
            [old[~old["tic_id"].isin(summary.get("tic_id", []))], summary], ignore_index=True
        )
    if len(summary):
        summary = summary.sort_values("tic_id").reset_index(drop=True)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        summary.to_parquet(out_path, index=False)
        okrows = summary[summary["status"] == "ok"]
        log.info(
            "wrote %s: %d ok, %d error; median runtime %.0f s/star; %d stars with SDE > 9 "
            "(%d of them with >= 2 transits with data)",
            out_path,
            len(okrows),
            (summary["status"] != "ok").sum(),
            okrows["runtime_s"].median(),
            (okrows["SDE"] > 9).sum(),
            ((okrows["SDE"] > 9) & (okrows["distinct_transit_count"] >= 2)).sum(),
        )
    return summary


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--sector", type=int, default=int(os.environ.get("SECTOR", 0)) or None
    )
    p.add_argument("--n-jobs", type=int, default=8, help="TLS is CPU-bound: use physical cores")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--force", action="store_true", help="re-search stars already in the output")
    p.add_argument(
        "--period-min", type=float, default=DEFAULT_PERIOD_MIN,
        help="shortest trial period in days (keep < 0.46 to recover LHS 3844 b)",
    )
    p.add_argument(
        "--period-max", type=float, default=None,
        help="longest trial period in days (default: time span / 2 per star)",
    )
    p.add_argument("--oversampling", type=int, default=3, help="TLS oversampling_factor (5 = careful)")
    p.add_argument(
        "--second-pass", action="store_true",
        help="mask each star's first-pass signal and search again; reads data/search/sector_N.parquet, "
        "writes data/search/pass2/sector_N_pass2.parquet with the same columns",
    )
    p.add_argument(
        "--only-candidates", action="store_true",
        help="with --second-pass: only stars whose first-pass signal passes vetting (sector_N_vetted.parquet)",
    )
    a = p.parse_args()
    if a.sector is None:
        p.error("no --sector given and SECTOR not set in .env")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    s = search_sector(
        a.sector, a.n_jobs, a.limit, a.force,
        a.period_min, a.period_max, a.oversampling, a.second_pass, a.only_candidates,
    )
    cols = [
        "tic_id", "status", "SDE", "snr", "period", "depth_ppm", "duration_h",
        "transit_count", "distinct_transit_count", "odd_even_mismatch", "runtime_s",
    ]
    print(s[[c for c in cols if c in s.columns]].to_string())


if __name__ == "__main__":
    main()
