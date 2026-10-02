"""Injection-recovery: how complete are steps 3-5 (detrend, search, vet)?

Usage:
    uv run python -m src.inject --n 1000 --n-jobs 11      # sector from .env; resumable
    uv run python -m src.inject --n 50 --seed 1           # a quick test

For each injection: pick a random star whose real signal did NOT pass vetting, read its raw light
curve with the same code as src.detrend (load_raw), multiply in a limb-darkened batman transit
of random period / planet radius / impact parameter / epoch, then run the real pipeline
functions on it: detrend.flatten_lc -> search.search_one -> vet.vet_one. Thresholds are NOT
applied here (they are applied in the analysis, so they can be tuned on the stored metrics).

Output: data/injection/sector_{N}_injrec.parquet — one row per injection: injected parameters,
expected SNR, the TLS result (search columns, snr -> snr_tls), every vet metric.
Written after every chunk, so an interrupted run resumes where it stopped.

Analysis and threshold tuning: notebooks/pipeline-8-injection-recovery.ipynb
"""

from __future__ import annotations

import argparse
import logging
import os
import time
from dataclasses import asdict
from pathlib import Path

import batman
import dotenv
import joblib
import numpy as np
import pandas as pd

from src.detrend import flatten_lc, load_raw
from src.search import DEFAULT_PERIOD_MIN, search_one
from src.vet import momentum_dumps, vet_one

log = logging.getLogger("inject")

dotenv.load_dotenv()
DATA_DIR = Path(
    os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data")
)
INJ_DIR = DATA_DIR / "injection"

G = 6.674e-11
R_SUN = 6.957e8
M_SUN = 1.989e30
DAY = 86400.0
R_EARTH_PER_R_SUN = 109.08
# Quadratic limb darkening for a ~3,400 K dwarf in the TESS band (Claret 2017 tables give
# u1 ~ 0.2-0.3, u2 ~ 0.3-0.4 for M dwarfs); the exact values barely matter for detectability.
LIMB_DARK = [0.25, 0.35]


def _worker_logging() -> None:
    """joblib (loky) workers are fresh processes with no logging config — set it up per worker."""
    if not logging.getLogger().handlers:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s [%(processName)s] %(message)s",
        )
        logging.getLogger("lightkurve").setLevel(logging.ERROR)
        logging.getLogger("search").setLevel(logging.WARNING)  # one line per injection is enough


def draw_injections(stars: pd.DataFrame, n: int, seed: int, p_range, rp_range) -> pd.DataFrame:
    """Random injection parameters. P and Rp log-uniform, b uniform in [0, 0.9], epoch uniform.

    Each injection has its own generator seeded by (seed, inj_id), so injection i is the same
    whatever n is — a run with --n 50 can be extended with --n 1000 and resumes consistently.
    """
    rows = []
    for i in range(n):
        rng = np.random.default_rng([seed, i])
        star = stars.iloc[rng.integers(0, len(stars))]
        rows.append({
            "inj_id": i,
            "tic_id": int(star["tic_id"]),
            "raw_file": star["raw_file"],
            "rad": float(star["rad"]),
            "mass": float(star["mass"]),
            "contratio": float(star["contratio"]),
            "P_inj": float(np.exp(rng.uniform(np.log(p_range[0]), np.log(p_range[1])))),
            "rp_inj_earth": float(np.exp(rng.uniform(np.log(rp_range[0]), np.log(rp_range[1])))),
            "b_inj": float(rng.uniform(0, 0.9)),
            "phase_inj": float(rng.uniform(0, 1)),  # epoch = first time + phase * P
        })
    return pd.DataFrame(rows)


def transit_model(t: np.ndarray, P: float, t0: float, rp_earth: float, b: float, rad: float, mass: float,
                  cadence_d: float) -> tuple[np.ndarray, float, float]:
    """batman flux model (1 out of transit); also returns the central-transit depth and duration [h]."""
    a = (G * mass * M_SUN * (P * DAY) ** 2 / (4 * np.pi**2)) ** (1 / 3) / (rad * R_SUN)  # a / R_star
    k = rp_earth / (rad * R_EARTH_PER_R_SUN)
    p = batman.TransitParams()
    p.t0, p.per, p.rp, p.a = t0, P, k, a
    p.inc = float(np.degrees(np.arccos(min(b / a, 1.0))))
    p.ecc, p.w = 0.0, 90.0
    p.limb_dark, p.u = "quadratic", LIMB_DARK
    # supersample over the exposure: 200-s QLP cadences smear a ~30-min transit noticeably
    flux = batman.TransitModel(p, t, supersample_factor=5, exp_time=cadence_d).light_curve(p)
    depth = 1 - float(batman.TransitModel(p, np.array([t0])).light_curve(p)[0])
    x = (1 + k) ** 2 - b**2
    dur_h = P / np.pi * np.arcsin(min(1.0, np.sqrt(max(x, 0)) / (a * np.sin(np.radians(p.inc))))) * 24
    return flux, depth, dur_h


def inject_one(inj: dict, dumps: np.ndarray, tmp_dir: Path) -> dict:
    """Inject one transit, run detrend -> search -> vet on it, return everything. Never raises."""
    _worker_logging()
    out = dict(inj)
    t_start = time.time()
    tmp = tmp_dir / f"{inj['inj_id']}.parquet"
    try:
        t, f, e, _ = load_raw(inj["raw_file"])
        cadence = float(np.median(np.diff(t)))
        t0 = t.min() + inj["phase_inj"] * inj["P_inj"]
        model, depth, dur_h = transit_model(
            t, inj["P_inj"], t0, inj["rp_inj_earth"], inj["b_inj"], inj["rad"], inj["mass"], cadence
        )
        flat, _ = flatten_lc(t, f * model, e)
        flat.to_parquet(tmp, index=False)

        # what an ideal matched filter would see: depth / point scatter * sqrt(points in transit)
        tt = flat["time"].to_numpy()
        dt = (tt - t0 + inj["P_inj"] / 2) % inj["P_inj"] - inj["P_inj"] / 2
        n_in = int((np.abs(dt) < dur_h / 48).sum())
        y = flat["flux"].to_numpy()
        sigma = 1.4826 * np.median(np.abs(y - np.median(y)))
        out.update(T0_inj=t0, depth_inj_ppm=depth * 1e6, dur_inj_h=dur_h, n_in_inj=n_in,
                   snr_expected=depth / sigma * np.sqrt(n_in) if sigma > 0 else np.nan)

        s = asdict(search_one(inj["tic_id"], str(tmp), inj["rad"], inj["mass"], DEFAULT_PERIOD_MIN))
        for k in ("transit_times", "transit_depths"):
            s.pop(k)
        s = {{"snr": "snr_tls", "status": "search_status", "error": "search_error",
              "runtime_s": "search_runtime_s", "n_points": "search_n_points"}.get(k, k): v for k, v in s.items()}
        out.update(s)
        if s["search_status"] == "ok":
            star = {"rad": inj["rad"], "mass": inj["mass"], "contratio": inj.get("contratio", np.nan)}
            v = asdict(vet_one(inj["tic_id"], s, star, str(tmp), inj["raw_file"], dumps))
            v.pop("tic_id")
            out.update(v)
        out["inj_runtime_s"] = time.time() - t_start
        log.info(
            "inj %d TIC %d: P %.3f d, Rp %.2f Re, depth %.0f ppm, SNR_exp %.1f -> TLS P %.3f d SDE %.1f",
            inj["inj_id"], inj["tic_id"], inj["P_inj"], inj["rp_inj_earth"], depth * 1e6,
            out["snr_expected"], out.get("period", np.nan), out.get("SDE", np.nan),
        )
    except Exception as ex:
        out["inj_error"] = f"{type(ex).__name__}: {ex}"
        out["inj_runtime_s"] = time.time() - t_start
        log.error("inj %d failed: %s", inj["inj_id"], out["inj_error"])
    finally:
        tmp.unlink(missing_ok=True)
    return out


def recovery_flags(df: pd.DataFrame, period_tol: float = 0.01) -> pd.DataFrame:
    """Add found (TLS found the injected period and epoch) and found_alias (P/2, 2P, ...)."""
    df = df.copy()
    ratio = df["period"] / df["P_inj"]
    dt = ((df["T0"] - df["T0_inj"]) / df["P_inj"] + 0.5) % 1 - 0.5
    on_time = (dt.abs() * df["P_inj"] * 24) < np.maximum(df["dur_inj_h"], 1.0)
    df["found"] = ((ratio - 1).abs() < period_tol) & on_time
    alias = np.zeros(len(df), bool)
    for r in (0.5, 2, 1 / 3, 3):
        alias |= ((ratio - r).abs() / r) < period_tol
    df["found_alias"] = alias & ~df["found"]
    return df


def inject_sector(
    sector: int,
    n: int = 1000,
    n_jobs: int = 8,
    seed: int = 0,
    p_range=(0.3, 12.0),
    rp_range=(0.5, 6.0),
    chunk: int = 55,
) -> pd.DataFrame:
    targets = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}.parquet",
                              columns=["tic_id", "rad", "mass", "contratio"])
    manifest = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}_manifest.parquet")
    detrend = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}_detrend.parquet")
    vetted = pd.read_parquet(DATA_DIR / "search" / f"sector_{sector}_vetted.parquet",
                             columns=["tic_id", "passes_all"])
    stars = (
        detrend[detrend["status"] == "ok"][["tic_id"]]
        .merge(manifest[manifest["status"] == "ok"][["tic_id", "file_path"]].rename(columns={"file_path": "raw_file"}),
               on="tic_id")
        .merge(targets, on="tic_id")
        .merge(vetted, on="tic_id")
    )
    # stars with a real vetted signal would confuse "recovered" with "found the real planet"
    stars = stars[~stars["passes_all"] & stars["rad"].notna() & stars["mass"].notna()].reset_index(drop=True)
    plan = draw_injections(stars, n, seed, p_range, rp_range)

    out_path = INJ_DIR / f"sector_{sector}_injrec.parquet"
    tmp_dir = INJ_DIR / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    done = pd.read_parquet(out_path) if out_path.exists() else pd.DataFrame(columns=["inj_id"])
    todo = plan[~plan["inj_id"].isin(done["inj_id"])]
    dumps = momentum_dumps(manifest)
    log.info("sector %d: %d stars eligible; %d injections planned, %d done, %d to run (%d workers)",
             sector, len(stars), n, len(done), len(todo), n_jobs)

    for start in range(0, len(todo), chunk):
        part = todo.iloc[start:start + chunk]
        rows = joblib.Parallel(n_jobs=n_jobs)(
            joblib.delayed(inject_one)(r._asdict(), dumps, tmp_dir) for r in part.itertuples(index=False)
        )
        done = pd.concat([done, pd.DataFrame(rows)], ignore_index=True) if len(done) else pd.DataFrame(rows)
        done = done.sort_values("inj_id").reset_index(drop=True)
        done.to_parquet(out_path, index=False)
        rec = recovery_flags(done[done.get("search_status", pd.Series(dtype=str)) == "ok"])
        log.info("progress %d/%d written; found so far %.0f %%", len(done), n, 100 * rec["found"].mean())
    return done


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--sector", type=int, default=int(os.environ.get("SECTOR", 0)) or None
    )
    p.add_argument("--n", type=int, default=1000, help="number of injections (resumes if the file exists)")
    p.add_argument("--n-jobs", type=int, default=8)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--p-min", type=float, default=0.3)
    p.add_argument("--p-max", type=float, default=12.0)
    p.add_argument("--rp-min", type=float, default=0.5, help="Earth radii")
    p.add_argument("--rp-max", type=float, default=6.0, help="Earth radii")
    a = p.parse_args()
    if a.sector is None:
        p.error("no --sector given and SECTOR not set in .env")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    logging.getLogger("lightkurve").setLevel(logging.ERROR)
    inject_sector(a.sector, a.n, a.n_jobs, a.seed, (a.p_min, a.p_max), (a.rp_min, a.rp_max))


if __name__ == "__main__":
    main()
