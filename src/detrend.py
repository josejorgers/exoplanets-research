"""Step 3 — detrend (flatten) the raw light curves with wotan's biweight filter.

Usage:
    uv run python -m src.detrend                      # sector from .env
    uv run python -m src.detrend --sector 104 --limit 20 --window-length 0.5
    uv run python -m src.detrend --edge-trim 0.25          # days dropped at each orbit edge (default)
    uv run python -m src.detrend --force

Input : data/targets/sector_{N}_manifest.parquet (rows with status == "ok", from src.fetch)
Output: data/lc_flat/{tic_id}.parquet  columns: time (BTJD), flux (flattened, ~1),
        flux_err, trend, flux_raw (normalized, before flattening)
        data/targets/sector_{N}_detrend.parquet   one summary row per star
"""

from __future__ import annotations

import argparse
import logging
import os
from dataclasses import asdict, dataclass
from pathlib import Path

import dotenv
import joblib
import lightkurve as lk
import numpy as np
import pandas as pd
import wotan

log = logging.getLogger("detrend")

dotenv.load_dotenv()
DATA_DIR = Path(
    os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data")
)
FLAT_DIR = DATA_DIR / "lc_flat"


@dataclass
class DetrendRow:
    tic_id: int
    n_in: int = 0  # points after quality mask + NaN removal
    n_out: int = 0  # points after upward outlier clipping
    n_trimmed: int = 0  # points dropped at the orbit edges (edge_trim)
    std_flat: float = np.nan  # robust std (1.4826*MAD) of flattened flux (ppm)
    med_err: float = np.nan  # median flux_err (ppm) — should be close to std_flat
    file_path: str = ""
    status: str = "error"  # ok | empty | error
    error: str = ""


def _worker_logging() -> None:
    """joblib (loky) workers are fresh processes with no logging config — set it up per worker."""
    if not logging.getLogger().handlers:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s [%(processName)s] %(message)s",
        )
        logging.getLogger("lightkurve").setLevel(logging.ERROR)


def edge_mask(time: np.ndarray, edge_trim: float, gap: float = 0.5) -> np.ndarray:
    """True for points farther than `edge_trim` days from the start/end of each orbit.

    Orbits are the segments separated by gaps longer than `gap` days (the mid-sector
    downlink). Scattered light as the spacecraft comes off perigee makes the raw flux
    climb ~10 % over the first day of an orbit; the biweight window is one-sided there
    and leaves a ~2 % residual dip that TLS reports as a high-SDE single "transit".
    Apply this AFTER flattening: trimming before only moves the edge (tested on
    TIC 6283622: 2.3 % residual untrimmed, 1.9 % trimmed before, 0.36 % = noise
    trimmed after). edge_trim = window_length / 2 covers the one-sided region.
    """
    keep = np.ones(len(time), dtype=bool)
    if edge_trim <= 0 or len(time) == 0:
        return keep
    starts = np.r_[0, np.where(np.diff(time) > gap)[0] + 1]
    ends = np.r_[starts[1:], len(time)]
    for a, b in zip(starts, ends):
        seg = time[a:b]
        keep[a:b] = (seg - seg[0] >= edge_trim) & (seg[-1] - seg >= edge_trim)
    return keep


def load_raw(raw_file: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Steps 1-2 of detrending: read, quality-mask, drop NaNs, normalize to 1, clip upward outliers.

    Returns (time [BTJD], flux, flux_err, n_in) where n_in counts the points before clipping.
    Also used by src.inject, which multiplies `flux` by a transit model before flatten_lc().
    """
    # 1. read (lightkurve applies the quality bitmask chosen at download), drop NaN flux, normalize to 1
    lc = lk.read(raw_file)
    # QLP flags bad cadences with bit 31 (2^30 = 1073741824), which lightkurve's default TESS
    # bitmask does not know about -> ~20 % junk points (flux ~0) survive unless we keep quality == 0.
    # SPOC files already had the "default" bitmask applied at download; leave them as they are.
    if str(lc.meta.get("AUTHOR", "")).upper() == "QLP":
        lc = lc[lc.quality.value == 0]
    lc = lc.remove_nans()
    n_in = len(lc)
    if n_in < 100:  # a few QLP files contain only NaN flux
        return np.array([]), np.array([]), np.array([]), n_in
    lc = lc.normalize()
    # 2. clip UPWARD outliers only (flares, cosmic rays). Never clip downward: those are the transits.
    lc = lc.remove_outliers(sigma_upper=4, sigma_lower=np.inf)
    # time in BTJD = BJD - 2457000
    return lc.time.value, lc.flux.value.astype(float), lc.flux_err.value.astype(float), n_in


def flatten_lc(
    time: np.ndarray,
    flux: np.ndarray,
    flux_err: np.ndarray,
    window_length: float = 0.5,
    edge_trim: float = 0.25,
) -> tuple[pd.DataFrame, int]:
    """Steps 3-4 of detrending: wotan biweight flatten, then trim orbit edges.

    Returns (DataFrame with time, flux, flux_err, trend, flux_raw; number of edge-trimmed points).
    """
    # 3. wotan returns plain numpy arrays: (flattened_flux, trend) with return_trend=True.
    flat, trend = wotan.flatten(
        time,
        flux,
        method="biweight",
        window_length=window_length,
        return_trend=True,
    )
    # biweight leaves NaN where the window has too few points (edges, gaps) — drop them
    good = np.isfinite(flat)
    # 4. drop the first/last edge_trim days of each orbit, where the window was one-sided (see edge_mask)
    keep = edge_mask(time, edge_trim)
    n_trimmed = int((good & ~keep).sum())
    good &= keep

    # QLP files carry NaN flux_err for the SAP flux; fall back to a robust scatter estimate
    bad_err = ~np.isfinite(flux_err)
    if bad_err.any():
        mad = np.nanmedian(np.abs(flat[good] - np.nanmedian(flat[good])))
        flux_err = np.where(bad_err, 1.4826 * mad, flux_err)

    out = pd.DataFrame(
        {
            "time": time[good],
            "flux": flat[good],
            "flux_err": flux_err[good],
            "trend": trend[good],
            "flux_raw": flux[good],
        }
    )
    return out, n_trimmed


def detrend_one(
    tic_id: int, raw_file: str, window_length: float = 0.5, edge_trim: float = 0.25
) -> DetrendRow:
    """Read one FITS light curve, normalize, clip flares, flatten with wotan, trim orbit edges, save parquet."""
    _worker_logging()
    row = DetrendRow(tic_id=tic_id)
    try:
        time, flux, flux_err, row.n_in = load_raw(raw_file)
        if len(time) == 0:
            row.status = "empty"
            log.warning("TIC %d: only %d finite points, skipping", tic_id, row.n_in)
            return row
        row.n_out = len(time)
        out, row.n_trimmed = flatten_lc(time, flux, flux_err, window_length, edge_trim)
        FLAT_DIR.mkdir(parents=True, exist_ok=True)
        path = FLAT_DIR / f"{tic_id}.parquet"
        out.to_parquet(path, index=False)

        # robust std (1.4826 * MAD) so a few leftover outliers do not dominate the metric
        row.std_flat = float(
            1.4826 * np.median(np.abs(out["flux"] - np.median(out["flux"]))) * 1e6
        )
        row.med_err = float(np.nanmedian(out["flux_err"]) * 1e6)
        row.file_path = str(path)
        row.status = "ok"
        log.info(
            "TIC %d: %d -> %d pts (%d edge-trimmed), std %.0f ppm, median err %.0f ppm",
            tic_id,
            row.n_in,
            len(out),
            row.n_trimmed,
            row.std_flat,
            row.med_err,
        )
    except Exception as e:
        row.error = f"{type(e).__name__}: {e}"
        log.error("TIC %d failed: %s", tic_id, row.error)
    return row


def detrend_sector(
    sector: int,
    window_length: float = 0.5,
    n_jobs: int = 8,
    limit: int | None = None,
    force: bool = False,
    edge_trim: float = 0.25,
) -> pd.DataFrame:
    manifest = pd.read_parquet(
        DATA_DIR / "targets" / f"sector_{sector}_manifest.parquet"
    )
    ok = manifest[manifest["status"] == "ok"].head(limit)
    if not force:
        ok = ok[~ok["tic_id"].map(lambda t: (FLAT_DIR / f"{t}.parquet").exists())]
    log.info(
        "sector %d: %d light curves downloaded, %d to detrend (window %.2f d, edge trim %.2f d, %d workers)",
        sector,
        (manifest["status"] == "ok").sum(),
        len(ok),
        window_length,
        edge_trim,
        n_jobs,
    )

    rows = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(detrend_one)(int(r.tic_id), r.file_path, window_length, edge_trim)
        for r in ok.itertuples()
    )
    summary = pd.DataFrame([asdict(r) for r in rows])
    summary_path = DATA_DIR / "targets" / f"sector_{sector}_detrend.parquet"
    if summary_path.exists():  # keep rows of stars not re-run (also with --force)
        old = pd.read_parquet(summary_path)
        summary = pd.concat(
            [old[~old["tic_id"].isin(summary.get("tic_id", []))], summary],
            ignore_index=True,
        )
    if len(summary):
        summary = summary.sort_values("tic_id").reset_index(drop=True)
        summary.to_parquet(summary_path, index=False)
        okrows = summary[summary["status"] == "ok"]
        log.info(
            "wrote %s: %d ok, %d error; median std/err ratio %.2f (expect ~1)",
            summary_path,
            len(okrows),
            (summary["status"] != "ok").sum(),
            (okrows["std_flat"] / okrows["med_err"]).median(),
        )
    return summary


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--sector", type=int, default=int(os.environ.get("SECTOR", 0)) or None
    )
    p.add_argument(
        "--window-length",
        type=float,
        default=0.5,
        help="biweight window in days (~3x longest transit)",
    )
    p.add_argument(
        "--edge-trim",
        type=float,
        default=0.25,
        help="days dropped at the start and end of each orbit after flattening (0 to disable; ~window/2)",
    )
    p.add_argument("--n-jobs", type=int, default=8)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    if a.sector is None:
        p.error("no --sector given and SECTOR not set in .env")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    logging.getLogger("lightkurve").setLevel(logging.ERROR)
    s = detrend_sector(a.sector, a.window_length, a.n_jobs, a.limit, a.force, a.edge_trim)
    print(s.to_string())


if __name__ == "__main__":
    main()
