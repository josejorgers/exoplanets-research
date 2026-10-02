"""Step 5 — automated vetting of the TLS signals: physical consistency tests, no early exit.

Usage:
    uv run python -m src.vet                           # sector from .env
    uv run python -m src.vet --sector 104 --limit 50 --n-jobs 8
    uv run python -m src.vet --prefilter               # skip light-curve tests for stars failing SDE or transit count
    uv run python -m src.vet --rethreshold --snr-min 6 # re-apply thresholds to the existing table, no recomputation
    uv run python -m src.vet --second-pass             # vet data/search/pass2/sector_{N}_pass2.parquet

Input : data/search/sector_{N}.parquet           (rows with status == "ok", from src.search)
        data/targets/sector_{N}.parquet          (rad, mass, contratio, from src.targets)
        data/targets/sector_{N}_detrend.parquet  (file_path of the flattened light curve)
        data/targets/sector_{N}_manifest.parquet (file_path of the raw FITS: centroids, momentum dumps)
Output: data/search/sector_{N}_vetted.parquet    one row per star: search columns + every metric,
        one pass_* per hard test, one flag_* per soft flag, passes_all, n_flags, failed, flags.
        Sorted by rank: passes_all first, then fewest flags, then highest snr.
        data/search/sector_{N}_vetted_meta.json  thresholds used and the sector's momentum dumps.

Workers compute only metrics; apply_thresholds() turns them into decisions, so thresholds can be
tuned with --rethreshold without touching a light curve. Convention: a test that cannot be
evaluated (metric NaN) passes — absence of evidence is not a rejection.

Walkthrough of every test, its physics and its threshold: notebooks/pipeline-5-vet-walkthrough.ipynb
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
import warnings
from dataclasses import asdict, dataclass
from pathlib import Path

import dotenv
import joblib
import numpy as np
import pandas as pd
import wotan
from astropy.io import fits
from astropy.timeseries import LombScargle
from scipy.optimize import OptimizeWarning, curve_fit
from transitleastsquares import transit_mask

log = logging.getLogger("vet")

dotenv.load_dotenv()
DATA_DIR = Path(
    os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data")
)
SEARCH_DIR = DATA_DIR / "search"

THRESH = dict(
    sde_min=9.0,  # TLS white-noise FAP ~1e-4
    snr_min=7.0,  # Kepler/TESS 7.1 sigma convention; TOI-122 b (snr 5.6 in s104) fails this
    n_transits_min=2,  # periodicity needs at least two events with data
    grid_edge_frac=0.98,  # period within 2 % of period_max = peak at the grid edge
    survive_frac=0.5,  # signal must keep >= 50 % of its sqrt(N)-expected SNR after dropping suspect transits
    odd_even_sigma_max=3.0,
    secondary_snr_max=5.0,  # 5 sigma: the scan tests many phases (look-elsewhere)
    dur_ratio_min=1 / 3,  # observed / expected central-transit duration
    dur_ratio_max=3.0,
    rp_earth_max=16.8,  # 1.5 R_Jup; giants around M dwarfs are not inflated
    var_harmonic_tol=0.005,  # period ratio within 0.5 % of an integer ...
    var_harmonic_n_max=10,
    var_amp_frac_max=0.15,  # ... and variability semi-amplitude >= 15 % of the depth
    v_shape_max=0.3,  # soft: trapezoid T23/T14 below this = V-shaped
    centroid_sigma_max=3.0,  # soft
    contratio_max=0.5,  # soft
)
# Used while computing metrics (not re-applicable with --rethreshold).
DUMP_WINDOW_D = 0.1  # +- days around a momentum dump (plus half a transit)

HARD = [
    "sde", "snr", "n_transits", "grid_edge", "single_event", "momentum_dump",
    "odd_even", "secondary", "duration", "radius", "variability",
]
SOFT = ["v_shaped", "centroid", "contamination"]

G = 6.674e-11
R_SUN = 6.957e8
M_SUN = 1.989e30
DAY = 86400.0
R_EARTH_PER_R_SUN = 109.08


@dataclass
class VetRow:
    tic_id: int
    n_points: int = 0
    # --- noise model ---
    sigma_ppm: float = np.nan  # robust point scatter, out of transit
    beta: float = np.nan  # red-noise factor on the transit-duration timescale (>= 1)
    # --- per-transit tests ---
    n_transits_used: int = 0  # transits with in-transit data
    snr: float = np.nan  # our beta-corrected combined SNR (inverse-variance weighted)
    max_ses_frac: float = np.nan  # strongest single transit SNR / combined SNR
    single_event_survival: float = np.nan  # SNR kept without the strongest transit, vs sqrt(N) expectation
    n_near_dump: int = 0
    dump_survival: float = np.nan  # same, without the transits near momentum dumps
    # --- odd/even, secondary ---
    depth_odd_ppm: float = np.nan
    depth_even_ppm: float = np.nan
    odd_even_sigma: float = np.nan
    secondary_snr: float = np.nan  # most significant dip away from the transit, any phase
    secondary_depth_ppm: float = np.nan
    secondary_phase: float = np.nan
    # --- the star ---
    duration_exp_h: float = np.nan  # central transit, circular orbit, Kepler's 3rd law
    duration_ratio: float = np.nan
    rp_earth: float = np.nan  # rp_rs * R_star, no dilution correction (lower bound if diluted)
    rp_earth_undiluted: float = np.nan  # * sqrt(1 + contratio)
    contratio: float = np.nan
    # --- variability ---
    var_period: float = np.nan  # Lomb-Scargle peak of the transit-masked flattened flux, days
    var_power: float = np.nan
    var_amp_ppm: float = np.nan
    var_period_ratio: float = np.nan  # P_tls / P_var
    var_amp_frac: float = np.nan  # variability semi-amplitude / transit depth
    # --- soft-flag metrics ---
    trap_depth_ppm: float = np.nan
    trap_t14_h: float = np.nan
    trap_flat_frac: float = np.nan  # T23/T14: 1 = box, 0 = V
    centroid_shift_px: float = np.nan
    centroid_sigma: float = np.nan
    # --- bookkeeping ---
    runtime_s: float = np.nan
    status: str = "error"  # ok | prefiltered | error
    error: str = ""  # star-level error, or "test: message; ..." for tests that failed to compute


def _worker_logging() -> None:
    """joblib (loky) workers are fresh processes with no logging config — set it up per worker."""
    if not logging.getLogger().handlers:
        logging.basicConfig(
            level=logging.INFO,
            format="%(asctime)s %(levelname)s [%(processName)s] %(message)s",
        )
    # empty-slice means and non-converging fits are handled per test (NaN); keep the log clean
    warnings.simplefilter("ignore", RuntimeWarning)
    warnings.simplefilter("ignore", OptimizeWarning)


def momentum_dumps(manifest: pd.DataFrame, n_files: int = 20, seed: int = 0) -> np.ndarray:
    """Times (BTJD) of reaction-wheel desaturations, from QUALITY bit 6 (value 32) of raw files.

    Spacecraft-wide, so the same for every star in the sector; SPOC and QLP use the same bit.
    """
    ok = manifest[manifest["status"] == "ok"]["file_path"]
    times = []
    for f in ok.sample(min(n_files, len(ok)), random_state=seed):
        try:
            d = fits.getdata(f, 1)
        except Exception as e:
            log.warning("cannot read %s for momentum dumps: %s", f, e)
            continue
        q = d["QUALITY"].astype(np.int64)
        times.append(d["TIME"][((q & 32) > 0) & np.isfinite(d["TIME"])])
    times = np.sort(np.concatenate(times)) if times else np.array([])
    if len(times) == 0:
        return times
    groups = np.split(times, np.where(np.diff(times) > 0.1)[0] + 1)  # one dump spans several cadences
    return np.array([float(g.mean()) for g in groups])


# ------------------------------------------------------------------------------------------------
# helpers: noise model, per-transit depths
# ------------------------------------------------------------------------------------------------


def _robust_std(x) -> float:
    """1.4826 * MAD: equals the standard deviation for Gaussian data, ignores outliers."""
    x = np.asarray(x)
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def _phase_time(t, period, T0):
    """Time from the nearest predicted mid-transit, days, in [-P/2, P/2)."""
    return (t - T0 + 0.5 * period) % period - 0.5 * period


def _noise_model(t, y, oot, dur) -> tuple[float, float]:
    """(sigma, beta): out-of-transit point scatter and the red-noise factor on timescale `dur`."""
    sigma = _robust_std(y[oot])
    tb, yb = t[oot], y[oot]
    edges = np.arange(tb.min(), tb.max() + dur, dur)
    idx = np.digitize(tb, edges)
    cnt = np.bincount(idx)
    tot = np.bincount(idx, yb)
    full = cnt >= max(3, 0.5 * np.median(cnt[cnt > 0]))  # skip bins cut by gaps
    if full.sum() < 5:
        return sigma, 1.0
    n_dur = float(np.median(cnt[full]))
    sigma_dur = _robust_std(tot[full] / cnt[full])
    return sigma, max(1.0, sigma_dur / (sigma / np.sqrt(n_dur)))


def _depth_err(y, mask, err1) -> tuple[float, float]:
    """Depth (1 - mean) of the points in mask and its uncertainty; err1 = beta * sigma."""
    n = int(mask.sum())
    if n == 0:
        return np.nan, np.nan
    return 1 - float(y[mask].mean()), err1 / np.sqrt(n)


def _per_transit(t, y, in_tr, epoch, T0, period, err1) -> pd.DataFrame:
    rows = []
    for k in np.unique(epoch[in_tr]):
        m = in_tr & (epoch == k)
        d, e = _depth_err(y, m, err1)
        rows.append((k, T0 + k * period, d, e))
    pt = pd.DataFrame(rows, columns=["epoch", "t_mid", "depth", "err"])
    pt["snr"] = pt["depth"] / pt["err"]
    return pt


def _combine(pt) -> float:
    """Inverse-variance weighted combined SNR of the transits in pt."""
    if len(pt) == 0:
        return np.nan
    w = 1 / pt["err"] ** 2
    return float((pt["depth"] * w).sum() / w.sum() * np.sqrt(w.sum()))


def _survival(pt, drop) -> float:
    """SNR kept after removing the rows in `drop`, relative to the sqrt(N) expectation (0 if < 2 kept)."""
    keep = pt[~drop]
    if len(keep) < THRESH["n_transits_min"]:
        return 0.0
    return float(_combine(keep) / (_combine(pt) * np.sqrt(len(keep) / len(pt))))


def _expected_duration_h(period_d, rad, mass, k=0.0) -> float:
    """Central-transit (b=0) duration for a circular orbit, hours. rad/mass in solar units."""
    a = (G * mass * M_SUN * (period_d * DAY) ** 2 / (4 * np.pi**2)) ** (1 / 3)
    return period_d / np.pi * np.arcsin(min(1.0, rad * R_SUN * (1 + k) / a)) * 24


def _trapezoid(x, x0, depth, t14, f):
    """x and durations in hours. f = T23/T14, 1 = box, 0 = V."""
    ax = np.abs(x - x0)
    t23 = f * t14
    y = np.ones_like(x)
    y[ax <= t23 / 2] = 1 - depth
    ingress = (ax > t23 / 2) & (ax < t14 / 2)
    y[ingress] = 1 - depth * (t14 / 2 - ax[ingress]) / ((t14 - t23) / 2 + 1e-12)
    return y


# ------------------------------------------------------------------------------------------------
# one star
# ------------------------------------------------------------------------------------------------


def vet_one(
    tic_id: int,
    sig: dict,
    star: dict,
    flat_file: str,
    raw_file: str | None,
    dumps: np.ndarray,
    prefilter: bool = False,
    mask: tuple[float, float, float] | None = None,
) -> VetRow:
    """Compute every vetting metric for one star's TLS signal. Never raises; decisions come later.

    sig  = search row (period, T0, duration_h, depth_ppm, rp_rs, SDE, distinct_transit_count)
    star = target row (rad, mass, contratio)
    mask = (period [d], duration [d], T0) of a first-pass signal to remove first (second pass).
    """
    _worker_logging()
    row = VetRow(tic_id=tic_id)
    t_start = time.time()
    errors = []

    def attempt(name, fn):
        try:
            fn()
        except Exception as e:
            errors.append(f"{name}: {type(e).__name__}: {e}")

    try:
        c = star.get("contratio", np.nan)
        row.contratio = float(c) if c is not None and np.isfinite(c) else np.nan
        P, T0, dur_h = float(sig["period"]), float(sig["T0"]), float(sig["duration_h"])
        dur = dur_h / 24
        depth = float(sig["depth_ppm"]) / 1e6

        # --- tests that use only the star (cheap, always computed) ---
        def star_tests():
            rad, mass = float(star["rad"]), float(star["mass"])
            if np.isfinite(rad) and np.isfinite(mass):
                row.duration_exp_h = _expected_duration_h(P, rad, mass, float(sig["rp_rs"]))
                row.duration_ratio = dur_h / row.duration_exp_h
                row.rp_earth = float(sig["rp_rs"]) * rad * R_EARTH_PER_R_SUN
                row.rp_earth_undiluted = row.rp_earth * np.sqrt(1 + np.nan_to_num(row.contratio))

        attempt("star", star_tests)

        if prefilter and (
            sig["SDE"] < THRESH["sde_min"] or sig["distinct_transit_count"] < THRESH["n_transits_min"]
        ):
            row.status = "prefiltered"
            row.runtime_s = time.time() - t_start
            return row

        # --- shared preparation: light curve, phases, masks, noise ---
        lc = pd.read_parquet(flat_file, columns=["time", "flux"])
        t, y = lc["time"].to_numpy(float), lc["flux"].to_numpy(float)
        good = np.isfinite(t) & np.isfinite(y)
        if mask is not None:
            good &= ~transit_mask(t, mask[0], 2 * mask[1], mask[2])
        t, y = t[good], y[good]
        row.n_points = len(t)
        dt = _phase_time(t, P, T0)
        in_tr = np.abs(dt) < dur / 2
        oot = np.abs(dt) > dur
        epoch = np.round((t - T0) / P).astype(int)
        sigma, beta = _noise_model(t, y, oot, dur)
        row.sigma_ppm, row.beta = sigma * 1e6, beta
        err1 = beta * sigma

        # --- per-transit: combined SNR, single-event dominance, momentum dumps ---
        def transit_tests():
            pt = _per_transit(t, y, in_tr, epoch, T0, P, err1)
            row.n_transits_used = len(pt)
            if len(pt) == 0:
                return
            row.snr = _combine(pt)
            row.max_ses_frac = float(pt["snr"].max() / row.snr)
            row.single_event_survival = _survival(pt, (pt.index == pt["snr"].idxmax()))
            window = DUMP_WINDOW_D + dur / 2
            near = (
                np.array([np.min(np.abs(tm - dumps)) < window for tm in pt["t_mid"]])
                if len(dumps) else np.zeros(len(pt), bool)
            )
            row.n_near_dump = int(near.sum())
            row.dump_survival = _survival(pt, near) if near.any() else 1.0

        attempt("transits", transit_tests)

        # --- odd/even ---
        def odd_even():
            d_o, e_o = _depth_err(y, in_tr & (epoch % 2 == 1), err1)
            d_e, e_e = _depth_err(y, in_tr & (epoch % 2 == 0), err1)
            row.depth_odd_ppm, row.depth_even_ppm = d_o * 1e6, d_e * 1e6
            row.odd_even_sigma = abs(d_o - d_e) / np.hypot(e_o, e_e)

        attempt("odd_even", odd_even)

        # --- secondary eclipse: most significant dip at any phase away from the transit ---
        def secondary():
            best = (0.0, np.nan, np.nan)
            for cph in np.arange(-P / 2, P / 2, dur / 4):
                if abs(cph) < 1.5 * dur or P - abs(cph) < 1.5 * dur:
                    continue
                m = np.abs((dt - cph + P / 2) % P - P / 2) < dur / 2
                if m.sum() < 3:
                    continue
                d, e = _depth_err(y, m, err1)
                if d / e > best[0]:
                    best = (d / e, d * 1e6, (cph / P) % 1)
            row.secondary_snr, row.secondary_depth_ppm, row.secondary_phase = best

        attempt("secondary", secondary)

        # --- variability: Lomb-Scargle of the transit-masked flattened flux ---
        def variability():
            ls = LombScargle(t[oot], y[oot])
            freq, power = ls.autopower(
                minimum_frequency=1 / 13, maximum_frequency=1 / 0.02, samples_per_peak=10
            )
            f_peak = freq[np.argmax(power)]
            model = ls.model(np.linspace(0, 1 / f_peak, 200), f_peak)
            amp = (model.max() - model.min()) / 2
            row.var_period = float(1 / f_peak)
            row.var_power = float(power.max())
            row.var_amp_ppm = float(amp * 1e6)
            row.var_period_ratio = float(P * f_peak)
            row.var_amp_frac = float(amp / depth) if depth > 0 else np.nan

        attempt("variability", variability)

        # --- shape: trapezoid fit to the fold ---
        def shape():
            x = dt * 24
            sel = np.abs(x) < 3 * dur_h
            popt, _ = curve_fit(
                _trapezoid, x[sel], y[sel], p0=[0, depth, dur_h, 0.5],
                bounds=([-dur_h / 2, 0, 0.3 * dur_h, 0], [dur_h / 2, 1, 3 * dur_h, 1]),
            )
            row.trap_depth_ppm, row.trap_t14_h, row.trap_flat_frac = popt[1] * 1e6, popt[2], popt[3]

        attempt("shape", shape)

        # --- centroid: flux-weighted centroid from the raw file, in vs out of transit ---
        def centroid():
            if not raw_file:
                return
            d = fits.getdata(raw_file, 1)
            names = d.columns.names
            cols = ("MOM_CENTR1", "MOM_CENTR2") if "MOM_CENTR1" in names else ("SAP_X", "SAP_Y")
            tt = np.asarray(d["TIME"], float)
            idx = np.clip(np.searchsorted(tt, t), 0, len(tt) - 1)
            same = np.abs(tt[idx] - t) < 1e-4  # match cadences by time (< 9 s)
            shift, err = [], []
            for col in cols:
                v = np.asarray(d[col], float)[idx]
                g = same & np.isfinite(v)
                if g.sum() < 10:  # some files carry an all-NaN centroid column
                    return
                trend = wotan.flatten(t[g], v[g], method="median", window_length=0.5, return_trend=True)[1]
                r = v[g] - trend
                fin = np.isfinite(r)
                it, oo = in_tr[g] & fin, oot[g] & fin
                if it.sum() == 0 or oo.sum() < 10:
                    return
                shift.append(r[it].mean() - np.median(r[oo]))
                err.append(_robust_std(r[oo]) / np.sqrt(it.sum()))
            row.centroid_shift_px = float(np.hypot(*shift))
            row.centroid_sigma = row.centroid_shift_px / float(np.hypot(*err) / np.sqrt(2))

        attempt("centroid", centroid)

        row.status = "ok"
        row.error = "; ".join(errors)
        row.runtime_s = time.time() - t_start
    except Exception as e:
        row.runtime_s = time.time() - t_start
        row.error = f"{type(e).__name__}: {e}"
        log.error("TIC %d failed: %s", tic_id, row.error)
    return row


# ------------------------------------------------------------------------------------------------
# decisions
# ------------------------------------------------------------------------------------------------


def apply_thresholds(df: pd.DataFrame, thresh: dict = THRESH) -> pd.DataFrame:
    """Add pass_* (hard), flag_* (soft), passes_all, n_flags, failed, flags; sort by rank.

    Works on the metric columns only, so it can be re-run on a saved table (--rethreshold).
    A NaN metric (test not evaluable, or not computed because prefiltered) passes and raises no flag.
    """
    df = df.copy()
    nan = lambda c: df[c].isna()  # noqa: E731
    df["pass_sde"] = df["SDE"] >= thresh["sde_min"]
    df["pass_snr"] = df["snr_tls"] >= thresh["snr_min"]
    df["pass_n_transits"] = df["distinct_transit_count"] >= thresh["n_transits_min"]
    df["pass_grid_edge"] = (df["period"] < thresh["grid_edge_frac"] * df["period_max"]) & np.isfinite(
        df["period_unc"]
    )
    df["pass_single_event"] = (df["single_event_survival"] >= thresh["survive_frac"]) | nan("single_event_survival")
    df["pass_momentum_dump"] = (df["dump_survival"] >= thresh["survive_frac"]) | nan("dump_survival")
    df["pass_odd_even"] = ~(df["odd_even_sigma"] > thresh["odd_even_sigma_max"])
    df["pass_secondary"] = ~(df["secondary_snr"] >= thresh["secondary_snr_max"])
    df["pass_duration"] = df["duration_ratio"].between(thresh["dur_ratio_min"], thresh["dur_ratio_max"]) | nan(
        "duration_ratio"
    )
    df["pass_radius"] = ~(df["rp_earth"] > thresh["rp_earth_max"])
    r = df["var_period_ratio"].where(df["var_period_ratio"] >= 1, 1 / df["var_period_ratio"])
    n = r.round()
    harmonic = (n >= 1) & (n <= thresh["var_harmonic_n_max"]) & ((r - n).abs() / r < thresh["var_harmonic_tol"])
    df["pass_variability"] = ~(harmonic & (df["var_amp_frac"] >= thresh["var_amp_frac_max"]))

    df["flag_v_shaped"] = df["trap_flat_frac"] < thresh["v_shape_max"]
    df["flag_centroid"] = df["centroid_sigma"] > thresh["centroid_sigma_max"]
    df["flag_contamination"] = df["contratio"] > thresh["contratio_max"]

    hard = df[[f"pass_{k}" for k in HARD]]
    soft = df[[f"flag_{k}" for k in SOFT]]
    errored = df["status"] == "error"
    df["passes_all"] = hard.all(axis=1) & ~errored
    df["n_flags"] = soft.sum(axis=1).astype(int)
    df["failed"] = hard.apply(lambda row: ",".join(k for k, ok in zip(HARD, row) if not ok), axis=1)
    df.loc[errored, "failed"] = "error"
    df["flags"] = soft.apply(lambda row: ",".join(k for k, f in zip(SOFT, row) if f), axis=1)
    return df.sort_values(["passes_all", "n_flags", "snr"], ascending=[False, True, False]).reset_index(drop=True)


def log_funnel(df: pd.DataFrame) -> None:
    """How many stars survive each hard test, applied in sequence."""
    alive = df["status"] != "error"
    log.info("funnel: %d stars vetted", alive.sum())
    for k in HARD:
        alive &= df[f"pass_{k}"]
        log.info("  ... and pass_%-14s %5d", k, alive.sum())


# ------------------------------------------------------------------------------------------------
# sector
# ------------------------------------------------------------------------------------------------


def output_path(sector: int, second_pass: bool = False) -> Path:
    if second_pass:
        return SEARCH_DIR / "pass2" / f"sector_{sector}_pass2_vetted.parquet"
    return SEARCH_DIR / f"sector_{sector}_vetted.parquet"


SEARCH_COLS = [
    "tic_id", "SDE", "snr", "period", "period_unc", "period_max", "T0", "duration_h",
    "depth_ppm", "rp_rs", "transit_count", "distinct_transit_count", "odd_even_mismatch",
]


def _write(df: pd.DataFrame, out_path: Path, dumps, thresh: dict) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, index=False)
    meta = {
        "thresholds": thresh,
        "dump_window_d": DUMP_WINDOW_D,
        "momentum_dumps_btjd": [float(x) for x in dumps],
        "hard_tests": HARD,
        "soft_flags": SOFT,
        "written": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    out_path.with_name(out_path.stem + "_meta.json").write_text(json.dumps(meta, indent=2))


def _log_survivors(df: pd.DataFrame, out_path: Path) -> None:
    okrows = df[df["status"] != "error"]
    log_funnel(df)
    log.info(
        "wrote %s: %d rows (%d ok, %d prefiltered, %d error, %d with test errors); median runtime %.2f s/star",
        out_path, len(df), (df["status"] == "ok").sum(), (df["status"] == "prefiltered").sum(),
        (df["status"] == "error").sum(), (df["error"].astype(bool) & (df["status"] == "ok")).sum(),
        okrows["runtime_s"].median(),
    )
    surv = df[df["passes_all"]]
    log.info("%d stars pass all hard tests", len(surv))
    for r in surv.itertuples():
        log.info(
            "  TIC %10d  P %8.4f d  depth %7.0f ppm  dur %5.2f h  Rp %5.1f Re  SNR %6.1f  SDE %5.1f  flags: %s",
            r.tic_id, r.period, r.depth_ppm, r.duration_h, r.rp_earth, r.snr, r.SDE, r.flags or "-",
        )


def vet_sector(
    sector: int,
    n_jobs: int = 8,
    limit: int | None = None,
    prefilter: bool = False,
    second_pass: bool = False,
    thresh: dict = THRESH,
) -> pd.DataFrame:
    search_path = (
        SEARCH_DIR / "pass2" / f"sector_{sector}_pass2.parquet" if second_pass
        else SEARCH_DIR / f"sector_{sector}.parquet"
    )
    search = pd.read_parquet(search_path)
    search = search[search["status"] == "ok"][SEARCH_COLS].rename(columns={"snr": "snr_tls"})
    targets = pd.read_parquet(
        DATA_DIR / "targets" / f"sector_{sector}.parquet", columns=["tic_id", "rad", "mass", "contratio"]
    ).set_index("tic_id")
    detrend = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}_detrend.parquet")
    flat = detrend[detrend["status"] == "ok"].set_index("tic_id")["file_path"]
    manifest = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}_manifest.parquet")
    raw = manifest[manifest["status"] == "ok"].set_index("tic_id")["file_path"]
    masks = {}
    if second_pass:
        # the flattened light curves still contain the first-pass transits: remove them as search.py did
        p1 = pd.read_parquet(SEARCH_DIR / f"sector_{sector}.parquet")
        p1 = p1[p1["status"] == "ok"].set_index("tic_id")
        masks = {
            tic: (float(r.period), float(r.duration_h) / 24, float(r.T0)) for tic, r in p1.iterrows()
        }

    stars = search[search["tic_id"].isin(flat.index)].head(limit)
    dumps = momentum_dumps(manifest)
    log.info(
        "sector %d%s: %d signals to vet (%d workers, prefilter %s); momentum dumps at BTJD %s",
        sector, " second pass" if second_pass else "", len(stars), n_jobs, prefilter,
        np.round(dumps, 3).tolist(),
    )

    empty_star = {"rad": np.nan, "mass": np.nan, "contratio": np.nan}
    rows = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(vet_one)(
            int(s.tic_id),
            s._asdict(),
            targets.loc[s.tic_id].to_dict() if s.tic_id in targets.index else empty_star,
            flat[s.tic_id],
            raw.get(s.tic_id),
            dumps,
            prefilter,
            masks.get(s.tic_id),
        )
        for s in stars.itertuples(index=False)
    )
    metrics = pd.DataFrame([asdict(r) for r in rows])
    df = stars.merge(metrics, on="tic_id", how="left")

    out_path = output_path(sector, second_pass)
    if out_path.exists():  # keep rows of stars not re-vetted in this run (e.g. after --limit)
        old = pd.read_parquet(out_path)
        df = pd.concat([old[~old["tic_id"].isin(df["tic_id"])], df], ignore_index=True)
    df = apply_thresholds(df[[c for c in df.columns if not c.startswith(("pass_", "flag_"))]], thresh)
    _write(df, out_path, dumps, thresh)
    _log_survivors(df, out_path)
    return df


def rethreshold(sector: int, second_pass: bool = False, thresh: dict = THRESH) -> pd.DataFrame:
    """Re-apply thresholds to an existing vetted table without recomputing any metric."""
    out_path = output_path(sector, second_pass)
    df = pd.read_parquet(out_path)
    meta_path = out_path.with_name(out_path.stem + "_meta.json")
    dumps = json.loads(meta_path.read_text())["momentum_dumps_btjd"] if meta_path.exists() else []
    drop = [c for c in df.columns if c.startswith(("pass_", "flag_"))] + ["passes_all", "n_flags", "failed", "flags"]
    df = apply_thresholds(df.drop(columns=drop), thresh)
    _write(df, out_path, dumps, thresh)
    _log_survivors(df, out_path)
    return df


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--sector", type=int, default=int(os.environ.get("SECTOR", 0)) or None
    )
    p.add_argument("--n-jobs", type=int, default=8)
    p.add_argument("--limit", type=int, default=None)
    p.add_argument(
        "--prefilter", action="store_true",
        help="skip the light-curve tests for stars that fail SDE or the transit count (they cannot pass)",
    )
    p.add_argument("--second-pass", action="store_true", help="vet the second-pass search table")
    p.add_argument(
        "--rethreshold", action="store_true",
        help="only re-apply thresholds to the existing vetted table (no light curves read)",
    )
    p.add_argument("--sde-min", type=float, default=THRESH["sde_min"])
    p.add_argument("--snr-min", type=float, default=THRESH["snr_min"])
    p.add_argument("--odd-even-max", type=float, default=THRESH["odd_even_sigma_max"])
    p.add_argument("--secondary-max", type=float, default=THRESH["secondary_snr_max"])
    a = p.parse_args()
    if a.sector is None:
        p.error("no --sector given and SECTOR not set in .env")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    thresh = {
        **THRESH,
        "sde_min": a.sde_min,
        "snr_min": a.snr_min,
        "odd_even_sigma_max": a.odd_even_max,
        "secondary_snr_max": a.secondary_max,
    }
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        if a.rethreshold:
            df = rethreshold(a.sector, a.second_pass, thresh)
        else:
            df = vet_sector(a.sector, a.n_jobs, a.limit, a.prefilter, a.second_pass, thresh)
    cols = [
        "tic_id", "passes_all", "n_flags", "SDE", "snr", "period", "depth_ppm", "duration_h",
        "rp_earth", "failed", "flags",
    ]
    print(df[df["passes_all"]][cols].to_string())


if __name__ == "__main__":
    main()
