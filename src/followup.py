"""Plan §12 — follow-up of individual candidates: is the signal real, on the target, and a planet?

Usage:
    uv run --extra followup python -m src.followup --tics 306101698,229804549
    uv run --extra followup python -m src.followup --tics 306101698 --skip-fit --skip-fpp   # quick checks only

For each TIC (ephemeris and star from data/search/sector_{N}_known.parquet and the target table):
  1. fetch_sectors       every TESS sector of the star (SPOC > TESS-SPOC > QLP per sector), processed
                         with the pipeline's own detrend.load_raw / flatten_lc
  2. sector_search       TLS in a narrow period window around P, sector by sector, and on all
                         sectors with and without the discovery sector: does the signal repeat?
  3. fit_transit         batman + emcee fit of the combined light curve (re-flattened with the
                         transits masked): Rp, b, a/R*, and the stellar density the transit implies
  4. difference_image    TESScut pixels: where on the sky does the flux drop during transit?
  5. gaia_neighbours     Gaia DR3 stars that are bright enough to cause the dip if they eclipsed
  6. fpp                 TRICERATOPS false-positive probability (Giacalone et al. 2021)
Output: data/followup/{tic}/ — multisector.parquet, sectors.csv, fit.json, fit_samples.npz,
        diffimage.png, gaia.csv, summary.json

Walkthrough with the theory of each check: notebooks/pipeline-9-candidate-followup.ipynb
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
import warnings
from pathlib import Path

import astropy.units as u
import batman
import dotenv
import numpy as np
import pandas as pd
import wotan
from astropy.coordinates import SkyCoord
from astropy.time import Time
from transitleastsquares import cleaned_array, transitleastsquares
from transitleastsquares.grid import period_grid

from src.background import trilegal_file
from src import neighbours as nb
from src.detrend import flatten_lc, load_raw

log = logging.getLogger("followup")

dotenv.load_dotenv()
DATA_DIR = Path(
    os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data")
)
FOLLOWUP_DIR = DATA_DIR / "followup"
PREFERENCE = {"SPOC": 0, "TESS-SPOC": 1, "QLP": 2}
G_CGS = 6.674e-8
R_SUN_CM, M_SUN_G = 6.957e10, 1.989e33
R_EARTH_PER_R_SUN = 109.08
PIXEL_ARCSEC = 21.0
TLS_MIN_GRID = 100  # transitleastsquares.tls_constants.MINIMUM_PERIOD_GRID_SIZE
# A sector can only refute a signal if the signal would have been seen there: expected SNR (discovery depth,
# that sector's scatter and number of in-transit points) at least this high, i.e. >= 2.5 sigma above the
# detection cut of repetition(). Weaker sectors (short baselines, 30-min cadence) are reported but not counted.
SNR_INFORMATIVE = 7.5


# ------------------------------------------------------------------------------------------------
# 1. all sectors
# ------------------------------------------------------------------------------------------------


def fetch_sectors(tic_id: int, out_dir: Path) -> pd.DataFrame:
    """Every sector of the star, one product per sector, detrended exactly like the pipeline."""
    import lightkurve as lk

    path = out_dir / "multisector.parquet"
    if path.exists():
        return pd.read_parquet(path)
    res = lk.search_lightcurve(f"TIC {tic_id}", author=list(PREFERENCE))
    tab = res.table.to_pandas()
    tab["i"] = np.arange(len(tab))
    tab["pref"] = tab["author"].map(PREFERENCE)
    best = tab.sort_values(["sequence_number", "pref", "exptime"]).groupby("sequence_number").head(1)
    parts = []
    for row in best.itertuples():
        try:
            lc = res[int(row.i)].download(download_dir=os.environ.get("LIGHTKURVE_CACHE_DIR"))
            t, f, e, _ = load_raw(lc.meta["FILENAME"])
            if len(t) == 0:
                continue
            flat, _ = flatten_lc(t, f, e)
            flat["sector"] = int(row.sequence_number)
            flat["author"] = row.author
            flat["cadence_s"] = float(row.exptime)
            parts.append(flat)
            log.info("TIC %d sector %d (%s, %.0f s): %d points", tic_id, row.sequence_number, row.author,
                     row.exptime, len(flat))
        except Exception as ex:  # a missing or broken file must not stop the other sectors
            log.warning("TIC %d sector %d: %s", tic_id, row.sequence_number, ex)
    out = pd.concat(parts, ignore_index=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_parquet(path, index=False)
    return out


# ------------------------------------------------------------------------------------------------
# 2. does it repeat?
# ------------------------------------------------------------------------------------------------


def window_oversampling(rad: float, mass: float, time_span: float, pmin: float, pmax: float, start: int = 3) -> int:
    """Smallest oversampling_factor >= start that puts >= TLS_MIN_GRID trial periods in [pmin, pmax].

    With fewer, TLS silently replaces the grid by its DEFAULT one (R = M = 1, no period limits), so a
    "narrow-window" search becomes a full search: one sector of a P = 4.5 d signal has only ~25 grid
    periods in +-2 %. The oversampling only refines the grid; the number of independent periods (and so
    the look-elsewhere effect) stays the same.
    """
    ovs = start
    while True:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # "very large grid" for multi-year baselines
            grid = period_grid(rad, mass, time_span, oversampling_factor=ovs)  # full grid: no fallback
        n = int(np.sum((grid > pmin) & (grid <= pmax)))
        if n >= TLS_MIN_GRID or ovs > 10_000:
            return ovs
        ovs = int(np.ceil(ovs * 1.2 * TLS_MIN_GRID / max(n, 1)))


def _tls_window(lc: pd.DataFrame, P: float, rad: float, mass: float, window: float, threads: int = 1):
    t, y, dy = cleaned_array(lc["time"].values, lc["flux"].values, lc["flux_err"].values)
    pmin, pmax = P * (1 - window), P * (1 + window)
    ovs = window_oversampling(rad, mass, float(t.max() - t.min()), pmin, pmax)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = transitleastsquares(t, y, dy).power(
            R_star=rad, M_star=mass, R_star_min=0.1, R_star_max=1.0, M_star_min=0.05, M_star_max=1.0,
            period_min=pmin, period_max=pmax, oversampling_factor=ovs,
            duration_grid_step=1.1, use_threads=threads, show_progress_bar=False, verbose=False,
        )
    if any("period_grid defaults" in str(w.message) for w in caught) or not (pmin <= r.period <= pmax):
        raise RuntimeError(f"TLS ignored the period window {pmin:.5f}-{pmax:.5f} d (returned {r.period:.5f})")
    return r, ovs


def expected_snr(lc: pd.DataFrame, P: float, dur_h: float, depth_ppm: float) -> float:
    """SNR a transit of depth_ppm would reach in lc: depth / point scatter * sqrt(in-transit points).

    In-transit points = all points * duration / P (no ephemeris needed: over a sector every phase is
    sampled about equally). Uses the discovery depth, which selection biases high (the detected SNR
    sits near the threshold partly by luck), so this is an upper estimate: a sector judged informative
    really is.
    """
    y = lc["flux"].to_numpy()
    sigma = 1.4826 * np.median(np.abs(y - np.median(y)))
    n_in = len(y) * dur_h / 24 / P
    return float(depth_ppm / 1e6 / sigma * np.sqrt(n_in)) if sigma > 0 else np.nan


def sector_search(ms: pd.DataFrame, P: float, rad: float, mass: float, discovery_sector: int,
                  window: float = 0.02, threads: int = 1, dur_h: float | None = None,
                  depth_ppm: float | None = None) -> pd.DataFrame:
    """TLS restricted to P*(1 +- window): per sector, all but the discovery sector, and all sectors.

    In a narrow window the look-elsewhere effect is small, so a real transit comes back at the same
    period with a significant SNR in the independent sectors; noise or a one-sector systematic does
    not. Judge by `snr`, not `SDE`: in a +-2 % window the "rest of the periodogram" that SDE compares
    against is mostly the wings of the peak itself, so SDE is not meaningful there.
    """
    groups = [(f"s{s:03d}", g) for s, g in ms.groupby("sector")]
    groups += [("all but discovery", ms[ms["sector"] != discovery_sector]), ("all sectors", ms)]
    rows = []
    for name, g in groups:
        if len(g) < 500:
            continue
        r, ovs = _tls_window(g, P, rad, mass, window, threads)
        snr_exp = expected_snr(g, P, dur_h, depth_ppm) if dur_h and depth_ppm else np.nan
        rows.append({
            "data": name, "n_points": len(g), "author": "/".join(sorted(g["author"].unique())),
            "oversampling": ovs, "snr_expected": snr_exp,
            "SDE": float(r.SDE), "period": float(r.period), "period_unc": float(r.period_uncertainty),
            "T0": float(r.T0), "depth_ppm": (1 - float(r.depth)) * 1e6, "snr": float(r.snr),
            "duration_h": float(r.duration) * 24, "n_transits": int(r.distinct_transit_count),
            "odd_even_mismatch": float(r.odd_even_mismatch),
        })
        log.info("%s: SDE %.1f P %.5f depth %.0f ppm (%d transits)", name, r.SDE, r.period, (1 - r.depth) * 1e6,
                 r.distinct_transit_count)
    return pd.DataFrame(rows)


def alias_search(ms: pd.DataFrame, P: float, rad: float, mass: float, window: float = 0.02,
                 threads: int = 1) -> pd.DataFrame:
    """Per sector, the same narrow search around 2P (rows "sNNN@2P").

    An eclipsing binary with two different eclipses, found at half its period, folds both
    eclipses together at P. At 2P only the deeper one is in the box, so the depth at 2P is clearly
    larger than at P. A planet has the same depth at P and 2P (every other transit). The ±2 %
    window at P can never return 2P itself, so this needs its own search.
    """
    rows = []
    for s, g in ms.groupby("sector"):
        if len(g) < 500 or (g["time"].max() - g["time"].min()) < 4 * P:  # need >= 2 transits at 2P
            continue
        try:
            r, ovs = _tls_window(g, 2 * P, rad, mass, window, threads)
        except Exception as ex:
            log.warning("sector %d at 2P: %s", s, ex)
            continue
        rows.append({"data": f"s{int(s):03d}@2P", "n_points": len(g), "author": "/".join(sorted(g["author"].unique())),
                     "oversampling": ovs, "SDE": float(r.SDE), "period": float(r.period),
                     "period_unc": float(r.period_uncertainty), "T0": float(r.T0),
                     "depth_ppm": (1 - float(r.depth)) * 1e6, "snr": float(r.snr),
                     "duration_h": float(r.duration) * 24, "n_transits": int(r.distinct_transit_count),
                     "odd_even_mismatch": float(r.odd_even_mismatch)})
    return pd.DataFrame(rows)


# A sector shows unequal eclipses when, detecting the signal at both P and 2P, the 2P depth is >= 1.2x the P depth
# AND the difference is >= 2 sigma (sigma = depth / snr for each). The best-fit depth at 2P is biased high at modest
# SNR (fewer transits, free phase), so the ratio alone fired on half the sectors of a signal that repeats at P.
ALIAS_DEPTH_RATIO = 1.2
ALIAS_SIGMA = 2.0


def repetition(sec: pd.DataFrame, P: float, depth_ppm: float, discovery_sector: int,
               tol: float = 0.005, snr_min: float = 5.0) -> dict:
    """Does the signal come back in the other sectors, one by one?

    TLS's SDE is meaningless in a narrow window, and its snr is unreliable on multi-sector data
    (negative values occur), so the test counts independent sectors whose OWN narrow search returns
    the same period (within single-sector precision, `tol`) with snr >= snr_min, and checks that their
    depths agree with the discovery depth. Sectors returning 2P or P/2 are counted separately: an
    eclipsing binary detected at half its period comes back at 2P where both eclipses are covered.
    That needs the "sNNN@2P" rows of alias_search(): a sector counts towards `repeats_at_alias` when
    it detects the signal at both P and 2P and the 2P depth is >= ALIAS_DEPTH_RATIO x the P depth, by
    >= ALIAS_SIGMA.

    Only INFORMATIVE sectors (snr_expected >= SNR_INFORMATIVE) can refute the signal: a sector where a
    real transit would reach SNR ~3 says nothing when it returns noise. A detection in any sector counts.
    When no single sector is informative, the narrow search on ALL OTHER SECTORS COMBINED decides, if it
    is informative itself: the signal repeats if it returns P (within `tol`) at a depth within a factor
    of 2 (TLS's snr is unreliable on multi-sector data, so it is not used there). `testable` is False
    when nothing is informative and nothing detected the signal: no verdict.
    A table without snr_expected (older runs) treats every sector as informative.
    """
    one = sec[sec["data"].str.match(r"s\d{3}$")].copy()
    at2 = sec[sec["data"].str.match(r"s\d{3}@2P$")].copy()
    one["sector"] = one["data"].str[1:].astype(int)
    other = one[one["sector"] != discovery_sector]
    snr_exp = other["snr_expected"] if "snr_expected" in other else pd.Series(np.inf, index=other.index)
    informative = other[snr_exp.fillna(np.inf) >= SNR_INFORMATIVE]
    same = other[((other["period"] - P).abs() / P < tol) & (other["snr"] >= snr_min)]
    deeper = []
    for _, r2 in at2.iterrows():
        r1 = other[other["data"] == r2["data"][:4]]
        if len(r1) and r1.index[0] in same.index and r2["snr"] >= snr_min and abs(r2["period"] - 2 * P) / (2 * P) < tol \
                and r2["depth_ppm"] >= ALIAS_DEPTH_RATIO * r1["depth_ppm"].iloc[0]:
            d1, s1 = r1["depth_ppm"].iloc[0], r1["snr"].iloc[0]
            diff_sigma = (r2["depth_ppm"] - d1) / np.hypot(d1 / s1, r2["depth_ppm"] / r2["snr"])
            if diff_sigma >= ALIAS_SIGMA:
                deeper.append(r2["data"])
    alias = pd.Series(deeper, dtype=object)
    ratio = float(np.median(same["depth_ppm"]) / depth_ppm) if len(same) else np.nan
    n_inf = len(informative)
    # half of the informative sectors, at least 2 when there are 2 or more (a single sector can be a fluke)
    needed = max(1, min(n_inf, max(2, int(np.ceil(n_inf / 2)))))
    per_sector = bool(len(same) >= needed and 0.5 < ratio < 2)
    comb = sec[sec["data"] == "all but discovery"]
    comb_inf, comb_rep, comb_ratio = False, False, np.nan
    if len(comb):
        c = comb.iloc[0]
        comb_inf = bool(c.get("snr_expected", np.inf) >= SNR_INFORMATIVE) if pd.notna(c.get("snr_expected", np.inf)) \
            else True
        comb_ratio = float(c["depth_ppm"] / depth_ppm)
        comb_rep = bool(abs(c["period"] - P) / P < tol and 0.5 < comb_ratio < 2)
    if n_inf > 0 or len(same) > 0:
        repeats, basis = per_sector, "sectors"
    elif comb_inf:
        repeats, basis = comb_rep, "combined"
    else:
        repeats, basis = False, "none"
    return {
        "n_other_sectors": len(other), "n_informative": n_inf, "n_same_period": len(same),
        "n_same_informative": int(informative.index.isin(same.index).sum()), "n_alias_period": len(alias),
        "depth_ratio": ratio, "combined_informative": comb_inf, "combined_repeats": comb_rep,
        "combined_depth_ratio": comb_ratio, "basis": basis,
        "testable": basis != "none",
        "repeats": repeats,
        "repeats_at_alias": bool(len(alias) >= 1 and len(alias) >= int(np.ceil(len(same) / 2))),
    }


# ------------------------------------------------------------------------------------------------
# 3. transit fit
# ------------------------------------------------------------------------------------------------


def masked_flatten(ms: pd.DataFrame, P: float, T0: float, dur_d: float, window_length: float = 0.5) -> pd.DataFrame:
    """Re-flatten each sector with the transits masked, so the trend cannot eat into their depth."""
    out = []
    for _, g in ms.groupby("sector"):
        t, fr = g["time"].to_numpy(), g["flux_raw"].to_numpy()
        dt = (t - T0 + P / 2) % P - P / 2
        oot = np.abs(dt) > 0.75 * dur_d
        trend = wotan.flatten(t[oot], fr[oot], method="biweight", window_length=window_length, return_trend=True)[1]
        good = np.isfinite(trend)
        g = g.copy()
        g["flux_fit"] = fr / np.interp(t, t[oot][good], trend[good])
        out.append(g)
    return pd.concat(out, ignore_index=True)


def _models(params_list, groups):
    """One batman model per cadence group (different supersampling)."""
    return {cad: batman.TransitModel(params_list, t, supersample_factor=max(1, int(round(cad / 40))),
                                     exp_time=cad / 86400) for cad, t in groups.items()}


def fit_transit(ms: pd.DataFrame, P: float, T0: float, dur_h: float, depth_ppm: float, rad: float, mass: float,
                nwalkers: int = 32, nsteps: int = 2500, seed: int = 0) -> tuple[dict, np.ndarray]:
    """batman + emcee fit. Parameters: T0, P, k = Rp/R*, b, log10(a/R*), q1, q2 (Kipping 2013 limb darkening).

    Uniform priors (k in [0, 0.5], b in [0, 1+k], a/R* in [1.5, 200], q1, q2 in [0, 1]), circular
    orbit, errors = per-sector robust scatter of the out-of-transit residuals.
    Returns a summary dict (median and 16/84 percentiles + derived quantities) and the flat chain.
    """
    import emcee

    dur = dur_h / 24
    d = masked_flatten(ms, P, T0, dur)
    dt = (d["time"] - T0 + P / 2) % P - P / 2
    d = d[np.abs(dt) < 3 * dur].copy()
    d["err"] = 0.0
    for s, g in d.groupby("sector"):  # per-sector white-noise level from the out-of-transit points
        gdt = (g["time"] - T0 + P / 2) % P - P / 2
        oot = g.loc[np.abs(gdt) > dur, "flux_fit"]
        d.loc[g.index, "err"] = 1.4826 * np.median(np.abs(oot - np.median(oot)))
    d = d[np.isfinite(d["flux_fit"]) & (d["err"] > 0)]
    # reference epoch near the middle of the data: decorrelates T0 and P
    n_mid = round((np.median(d["time"]) - T0) / P)
    T0_mid = T0 + n_mid * P
    groups = {cad: g["time"].to_numpy() for cad, g in d.groupby("cadence_s")}
    obs = {cad: (g["flux_fit"].to_numpy(), g["err"].to_numpy()) for cad, g in d.groupby("cadence_s")}

    pm = batman.TransitParams()
    pm.t0, pm.per, pm.rp, pm.a, pm.inc, pm.ecc, pm.w = T0_mid, P, 0.05, 20, 89, 0, 90
    pm.limb_dark, pm.u = "quadratic", [0.3, 0.3]
    models = _models(pm, groups)

    def lnprob(theta):
        t0, per, k, b, log_a, q1, q2 = theta
        a = 10**log_a
        if not (0 < k < 0.5 and 0 <= b < 1 + k and 1.5 < a < 200 and 0 < q1 < 1 and 0 < q2 < 1):
            return -np.inf
        if abs(per - P) > 0.01 * P or abs(t0 - T0_mid) > 0.5 * dur:
            return -np.inf
        pm.t0, pm.per, pm.rp, pm.a = t0, per, k, a
        pm.inc = float(np.degrees(np.arccos(b / a)))
        sq = np.sqrt(q1)
        pm.u = [2 * sq * q2, sq * (1 - 2 * q2)]
        chi2 = 0.0
        for cad, m in models.items():
            f, e = obs[cad]
            chi2 += np.sum(((f - m.light_curve(pm)) / e) ** 2)
        return -0.5 * chi2

    a0 = (6.674e-11 * mass * 1.989e30 * (P * 86400) ** 2 / (4 * np.pi**2)) ** (1 / 3) / (rad * 6.957e8)
    p0 = np.array([T0_mid, P, np.sqrt(max(depth_ppm, 100) / 1e6), 0.3, np.log10(a0), 0.4, 0.3])
    rng = np.random.default_rng(seed)
    scale = np.array([dur / 20, P * 1e-5, 0.003, 0.1, 0.03, 0.1, 0.1])
    start = p0 + scale * rng.normal(size=(nwalkers, len(p0)))
    start[:, 3] = np.abs(start[:, 3])
    start[:, 5:] = np.clip(start[:, 5:], 0.05, 0.95)
    sampler = emcee.EnsembleSampler(nwalkers, len(p0), lnprob)
    t_start = time.time()
    sampler.run_mcmc(start, nsteps, progress=False)
    chain = sampler.get_chain(discard=nsteps // 3, thin=5, flat=True)
    names = ["T0", "P", "k", "b", "log10_a_rs", "q1", "q2"]
    summ = {"n_points": int(len(d)), "sectors": sorted(int(s) for s in d["sector"].unique()),
            "acceptance": float(np.mean(sampler.acceptance_fraction)), "runtime_s": time.time() - t_start}
    for i, nme in enumerate(names):
        lo, med, hi = np.percentile(chain[:, i], [16, 50, 84])
        summ[nme] = {"median": float(med), "minus": float(med - lo), "plus": float(hi - med)}
    # derived: planet radius (TIC radius with an assumed 7 % uncertainty), stellar density, duration
    rad_draw = rad * (1 + 0.07 * rng.normal(size=len(chain)))
    rp = chain[:, 2] * rad_draw * R_EARTH_PER_R_SUN
    a_rs = 10 ** chain[:, 4]
    rho_transit = 3 * np.pi * a_rs**3 / (G_CGS * (chain[:, 1] * 86400) ** 2)  # g/cm^3
    rho_tic = mass * M_SUN_G / (4 / 3 * np.pi * (rad * R_SUN_CM) ** 3)
    for nme, arr in [("rp_earth", rp), ("a_rs", a_rs), ("rho_star_transit", rho_transit)]:
        lo, med, hi = np.percentile(arr, [16, 50, 84])
        summ[nme] = {"median": float(med), "minus": float(med - lo), "plus": float(hi - med)}
    summ["rho_star_tic"] = float(rho_tic)
    summ["rho_ratio_median"] = float(np.median(rho_transit) / rho_tic)
    return summ, chain


# ------------------------------------------------------------------------------------------------
# 4. pixels: difference image
# ------------------------------------------------------------------------------------------------


def target_coord(star: pd.Series, epoch_btjd: float) -> SkyCoord:
    """TIC position (J2000) moved by the proper motion to the observation epoch."""
    years = (Time(epoch_btjd + 2457000, format="jd").jyear - 2000.0)
    pmra = (star.get("pmRA", 0) or 0) if np.isfinite(star.get("pmRA", np.nan)) else 0.0
    pmdec = (star.get("pmDEC", 0) or 0) if np.isfinite(star.get("pmDEC", np.nan)) else 0.0
    dec = star.dec + pmdec / 3.6e6 * years
    ra = star.ra + pmra / 3.6e6 * years / np.cos(np.radians(star.dec))
    return SkyCoord(ra, dec, unit="deg")


def _event_images(t: np.ndarray, cube: np.ndarray, centres, dur: float, exclude: np.ndarray | None = None):
    """Per event: median out-of-transit image minus mean in-transit image, around each time in `centres`.

    In-transit is the whole transit (|dt| < dur/2), out-of-transit the 0.75-2.5 durations on both sides.
    Cadences in `exclude` (the real transits, when building pseudo-events) are never used.
    """
    diffs, oots = [], []
    for tc in centres:
        dt = np.abs(t - tc)
        i_in, i_oot = dt < 0.5 * dur, (dt > 0.75 * dur) & (dt < 2.5 * dur)
        if exclude is not None:
            i_in, i_oot = i_in & ~exclude, i_oot & ~exclude
        if i_in.sum() >= 3 and i_oot.sum() >= 6:
            oot = np.median(cube[i_oot], axis=0)
            diffs.append(oot - np.mean(cube[i_in], axis=0))
            oots.append(oot)
    return np.array(diffs), np.array(oots)


DI_N_NULL = 60  # pseudo-transit difference images for the noise estimate
DI_N_NULL_MIN = 15  # fewer than this (P < ~8 durations): fall back to the transit-to-transit scatter
DI_SNR_MIN = 5.0  # a flux loss anywhere in the cutout counts only above max(this, the null's 99th percentile)
DI_SNR_MIN_TARGET = 3.0  # same, for the loss within 3 px of the target


def difference_image_stats(t: np.ndarray, cube: np.ndarray, P: float, T0: float, dur: float, tx: float, ty: float,
                           n_boot: int = 200, seed: int = 0) -> dict:
    """Difference image of a background-subtracted cube (time, y, x) and where its flux loss sits.

    Noise: the transit-to-transit scatter of each pixel is a poor estimate with few transits (with 4-6
    transits a few of 225 pixels get a tiny scatter by chance and a S/N of 10+; TIC 372868160). Instead,
    build the same difference image at DI_N_NULL "pseudo-transits" (the same ephemeris shifted to times with no
    transit, real transits masked). Their pixel-by-pixel spread is the noise of the real image, including
    systematics; their largest S/N is what pure noise produces in this cutout, so a real loss must beat it.
    The centroid error is the larger of the bootstrap over transits and the centroid spread when the
    pseudo-images are added to the real one as noise.
    """
    epoch = np.round((t - T0) / P)
    centres = T0 + np.unique(epoch) * P
    diffs, oots = _event_images(t, cube, centres, dur)
    if len(diffs) < 2:
        return {"n_transits": len(diffs), "error": "fewer than 2 usable transits in the cutout"}
    diff, oot = diffs.mean(axis=0), oots.mean(axis=0)
    yy, xx = np.indices(diff.shape)
    near_target = np.hypot(xx - tx, yy - ty) <= 3.0

    def centroid(img):
        """Flux-weighted centroid of the positive pixels within 1.5 px of the image's peak inside 3 px
        of the target. Restricting to the target's surroundings keeps noise peaks elsewhere in the
        cutout (and brighter neighbours in the out-of-transit image) from hijacking the centroid."""
        masked = np.where(near_target, img, -np.inf)
        py, px = np.unravel_index(np.argmax(masked), img.shape)
        w = np.where((np.hypot(xx - px, yy - py) <= 1.5) & (img > 0) & near_target, img, 0.0)
        if w.sum() <= 0:
            return np.nan, np.nan
        return float((w * xx).sum() / w.sum()), float((w * yy).sum() / w.sum())

    # pseudo-transits: shifts of 3 durations to P - 3 durations, avoiding phase 0.5 (an EB's secondary)
    real = np.abs(t - (T0 + epoch * P)) < 0.75 * dur
    shifts = np.linspace(3 * dur, P - 3 * dur, DI_N_NULL) if P > 6 * dur else np.array([])
    shifts = shifts[np.abs(shifts - P / 2) > 1.5 * dur]
    null = []
    for s in shifts:
        d_k, _ = _event_images(t, cube, centres + s, dur, exclude=real)
        if len(d_k) >= max(2, len(diffs) // 2):
            null.append(d_k.mean(axis=0) * np.sqrt(len(d_k) / len(diffs)))  # scaled to the real number of events
    null = np.array(null)
    if len(null) >= DI_N_NULL_MIN:
        noise_src = "pseudo-transits"
        sigma = null.std(axis=0, ddof=1)
        null_snr = null / sigma
        null_max = float(np.quantile(null_snr.reshape(len(null), -1).max(axis=1), 0.99))
        null_max_near = float(np.quantile(np.where(near_target, null_snr, -np.inf).reshape(len(null), -1).max(axis=1),
                                          0.99))
    else:
        noise_src = "transit-to-transit scatter"
        sigma = diffs.std(axis=0, ddof=1) / np.sqrt(len(diffs))
        null_max, null_max_near = np.nan, np.nan
    snr_img = diff / sigma

    cx, cy = centroid(diff)
    ox, oy = centroid(oot)
    rng = np.random.default_rng(seed)
    boot = np.array([centroid(diffs[rng.integers(0, len(diffs), len(diffs))].mean(axis=0)) for _ in range(n_boot)])
    boot = boot[np.all(np.isfinite(boot), axis=1)]
    sx, sy = boot.std(axis=0) if len(boot) > 10 else (np.nan, np.nan)
    if len(null) >= DI_N_NULL_MIN:
        pn = np.array([centroid(diff + n) for n in null])
        pn = pn[np.all(np.isfinite(pn), axis=1)]
        if len(pn) > 10:
            nx, ny = pn.std(axis=0)
            sx, sy = np.fmax(sx, nx), np.fmax(sy, ny)
    sig = float(np.hypot(sx, sy))
    off = float(np.hypot(cx - ox, cy - oy))  # difference centroid vs the target's own light centroid
    off_cat = float(np.hypot(cx - tx, cy - ty))  # vs the catalog position
    gy, gx = np.unravel_index(np.nanargmax(snr_img), snr_img.shape)
    ix, iy = int(round(float(tx))), int(round(float(ty)))
    ny_, nx_ = diff.shape
    near_snr = float(np.nanmax(np.where(near_target, snr_img, -np.inf)))
    thr_any = max(DI_SNR_MIN, null_max) if np.isfinite(null_max) else DI_SNR_MIN
    thr_near = max(DI_SNR_MIN_TARGET, null_max_near) if np.isfinite(null_max_near) else DI_SNR_MIN_TARGET
    res = {
        "n_transits": len(diffs), "target_px": [float(tx), float(ty)],
        "diff_centroid_px": [cx, cy], "oot_centroid_px": [ox, oy],
        "offset_px": off, "offset_arcsec": off * PIXEL_ARCSEC, "offset_err_px": sig,
        "offset_sigma": off / sig if sig > 0 else np.nan, "offset_from_catalog_px": off_cat,
        "strongest_loss_px": [int(gx), int(gy)], "strongest_loss_snr": float(snr_img[gy, gx]),
        "strongest_loss_dist_px": float(np.hypot(gx - tx, gy - ty)),
        "loss_snr_at_target": float(snr_img[iy, ix]) if 0 <= ix < nx_ and 0 <= iy < ny_ else np.nan,
        "loss_snr_near_target": near_snr,
        "noise_source": noise_src, "n_null": int(len(null)),
        "null_max_snr": null_max, "null_max_snr_near": null_max_near,
        "snr_threshold": thr_any, "snr_threshold_near": thr_near,
    }
    detected_near = near_snr >= thr_near
    detected_any = res["strongest_loss_snr"] >= thr_any
    # off target: the centroid of the flux loss is significantly away from the star's own light AND from its
    # catalog position (and there is a significant loss near the target to take a centroid of), or the
    # strongest loss in the cutout sits > 1.5 px away, beats what noise gives, and is clearly stronger than at
    # the target. Both positions are needed: in a crowded field the out-of-transit light centroid is pulled
    # towards the neighbours (TIC 214311908: 0.76 px from the light centroid, 0.30 px from the catalog position).
    centroid_off = detected_near and res["offset_sigma"] > 3 and min(off, off_cat) > 0.5
    peak_off = (detected_any and res["strongest_loss_dist_px"] > 1.5
                and res["loss_snr_at_target"] < 0.5 * res["strongest_loss_snr"])
    res["off_target"] = bool(centroid_off or peak_off)
    res["status"] = ("off_target" if res["off_target"] else
                     "on_target" if detected_near or detected_any else "inconclusive")
    res["_arrays"] = dict(diff=diff, oot=oot, target=[tx, ty], centroid=[cx, cy], oot_centroid=[ox, oy],
                          strongest=[gx, gy], snr_img=snr_img, boot=boot, sigma=sigma)
    return res


def difference_image(tic_id: int, star: pd.Series, P: float, T0: float, dur_h: float, sector: int,
                     out_dir: Path, size: int = 15, n_boot: int = 200, seed: int = 0) -> dict:
    """Mean out-of-transit minus mean in-transit image, transit by transit, from a TESScut cutout.

    The difference image shows where the flux disappears during transit. Its centroid is compared
    with the target's (proper-motion corrected) position; see difference_image_stats for the noise model
    and the on target / off target / inconclusive decision.
    """
    import lightkurve as lk

    tpf = lk.search_tesscut(f"TIC {tic_id}", sector=sector)[0].download(
        cutout_size=size, download_dir=os.environ.get("LIGHTKURVE_CACHE_DIR"))
    good = (tpf.quality == 0) & np.all(np.isfinite(tpf.flux.value), axis=(1, 2))
    t, cube = tpf.time.value[good], tpf.flux.value[good]
    # Per-cadence sky background: scattered light changes by more than a transit within an hour and
    # differently from one transit to the next; without this the difference image is dominated by it.
    # The 20th percentile of the cutout is a robust sky estimate even in crowded fields.
    cube = cube - np.percentile(cube.reshape(len(cube), -1), 20, axis=1)[:, None, None]
    tgt = target_coord(star, float(np.median(t)))
    tx, ty = tpf.wcs.world_to_pixel(tgt)
    res = {"sector": sector, **difference_image_stats(t, cube, P, T0, dur_h / 24, float(tx), float(ty),
                                                      n_boot=n_boot, seed=seed)}
    if "error" in res:
        return res
    gx, gy = res["strongest_loss_px"]
    sky = tpf.wcs.pixel_to_world(gx, gy)
    res["strongest_loss_radec"] = [float(sky.ra.deg), float(sky.dec.deg)]
    np.savez_compressed(out_dir / f"diffimage_s{sector}.npz", **res.pop("_arrays"))
    res["_tpf_wcs"] = tpf.wcs  # for plotting in the notebook; dropped from JSON
    return res


# ------------------------------------------------------------------------------------------------
# 5. Gaia neighbours
# ------------------------------------------------------------------------------------------------


def gaia_neighbours(coord: SkyCoord, depth_ppm: float, radius_arcsec: float = 63.0) -> pd.DataFrame:
    """Gaia DR3 stars within `radius` (default 3 TESS pixels) and whether they could cause the dip.

    A neighbour fainter by dm can produce at most a dip of 10^(-0.4 dm) of the total light (if it
    eclipsed completely). So only stars with dm < -2.5 log10(depth) are possible sources.
    """
    from astroquery.vizier import Vizier

    res = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Gmag", "BP-RP", "Plx"], row_limit=-1).query_region(
        coord, radius=radius_arcsec * u.arcsec, catalog="I/355/gaiadr3")
    if not len(res):
        return pd.DataFrame()
    g = res[0].to_pandas()
    g["sep_arcsec"] = SkyCoord(g["RA_ICRS"].values, g["DE_ICRS"].values, unit="deg").separation(coord).arcsec
    g = g.sort_values("sep_arcsec").reset_index(drop=True)
    g["dG"] = g["Gmag"] - g["Gmag"].iloc[0]  # relative to the nearest (= target) source
    dm_max = -2.5 * np.log10(depth_ppm / 1e6)
    g["could_cause_dip"] = g["dG"] < dm_max
    return g


# ------------------------------------------------------------------------------------------------
# 6. false-positive probability
# ------------------------------------------------------------------------------------------------


def fpp(tic_id: int, ms: pd.DataFrame, P: float, T0: float, dur_h: float, depth_ppm: float, sector: int,
        n_draws: int = 20_000, seed: int = 0, n_runs: int = 5, background: bool = True,
        exclusions: bool = True) -> dict:
    """TRICERATOPS false-positive probability from the discovery sector's folded light curve.

    TRICERATOPS integrates each scenario by Monte Carlo, and the scenario probabilities move a lot from
    run to run (for TIC 214311908 one neighbour's NTP went from 0.04 to 0.59 between two seeds), so
    calc_probs is repeated n_runs times and the mean and scatter are reported.

    Neighbours (src/neighbours.py): TIC DUPLICATE entries get no light and no scenarios; neighbours that
    could cause the dip but lack TIC Teff/radius/mass get Gaia DR3 values instead of TRICERATOPS' solar
    defaults; with exclusions=True the stars in data/followup/<tic>/pixel_excluded.json (ruled out by the
    per-pixel PRF fit) keep their light but get no scenarios."""
    from src import _compat  # noqa: F401  (scipy/pkg_resources shims needed by pytransit)
    import triceratops.triceratops as tr

    g = ms[ms["sector"] == sector]
    dt = (g["time"].to_numpy() - T0 + P / 2) % P - P / 2
    sel = np.abs(dt) < 3 * dur_h / 24
    oot = np.abs(dt) > dur_h / 24
    # Background population for the BTP/BEB/DTP/DEB scenarios: our own TRILEGAL query (see
    # src/background.py for why triceratops' built-in query always failed), cached per target.
    star = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}.parquet").set_index("tic_id").loc[tic_id]
    try:
        if not background:
            raise RuntimeError("background scenarios switched off by the caller")
        trilegal = trilegal_file(tic_id, float(star["ra"]), float(star["dec"]), FOLLOWUP_DIR / str(tic_id))
        trilegal_error = None
    except Exception as ex:
        trilegal, trilegal_error = None, f"{type(ex).__name__}: {ex}"
        log.warning("TIC %d: TRILEGAL query failed (%s)", tic_id, trilegal_error)
    # with trilegal_fname=None triceratops tries its own (broken) query; a missing file string stops it
    target = tr.target(ID=tic_id, sectors=np.array([sector]), trilegal_fname=trilegal or "none")
    # Aperture: pixels within 1.5 px of the target, close to QLP/SPOC apertures for T ~ 11-14
    # (triceratops' default 5x5 aperture lets far more neighbour light in than our light curves had).
    tx, ty = np.round(target.pix_coords[0][0])
    xx, yy = np.meshgrid(np.arange(tx - 2, tx + 3), np.arange(ty - 2, ty + 3))
    keep = np.hypot(xx - target.pix_coords[0][0][0], yy - target.pix_coords[0][0][1]) <= 1.5
    ap = np.column_stack([xx[keep], yy[keep]])
    # NB: despite its docstring ("[ppm]"), calc_depths expects a FRACTION: it computes depth / flux_ratio
    # and discards anything > 1, so passing ppm silently removes every star (all scenarios get prob 1/12).
    # Our light curves are DEBLENDED (neighbour light removed); triceratops models the light curve as
    # observed in the aperture, diluted by every star in it. So we first get the target's flux ratio in
    # the aperture, then re-dilute the light curve and the depth by it.
    st = target.stars
    dup = [i for i in st.index[1:] if st.at[i, "disposition"] == "DUPLICATE"]
    duplicates = [str(st.at[i, "ID"]) for i in dup if st.at[i, "Tmag"] < 20]
    st.loc[dup, "Tmag"] = nb.NO_FLUX_TMAG  # a second copy of another TIC star: no light of its own
    target.calc_depths(tdepth=0.5 * depth_ppm / 1e6, all_ap_pixels=[ap])
    fr_target = float(target.stars["fluxratio"].iloc[0])
    target.calc_depths(tdepth=depth_ppm / 1e6 * fr_target, all_ap_pixels=[ap])
    if not (target.stars["tdepth"] > 0).any():
        raise RuntimeError("no star can produce the transit depth; check the depth units")
    st = target.stars
    rows = [i for i in st.index[1:] if st.at[i, "tdepth"] > 0]  # neighbours that get NTP/NEB scenarios
    nb_log, nb_error, gaia_ids = [], None, pd.Series(pd.NA, index=rows, dtype="Int64")
    if rows:
        try:
            gaia = nb.query_gaia(float(st["ra"].iloc[0]), float(st["dec"].iloc[0]),
                                 float(st.loc[rows, "sep (arcsec)"].max()) + 5.0)
            gaia_ids = nb.match_gaia(st.loc[rows], gaia)
            nb_log = nb.fill_parameters(st, rows, gaia)
        except Exception as ex:
            nb_error = f"{type(ex).__name__}: {ex}"
            log.warning("TIC %d: Gaia parameters for neighbours failed (%s); solar defaults kept", tic_id, nb_error)
    excl = nb.load_exclusions(FOLLOWUP_DIR / str(tic_id) / "pixel_excluded.json") if exclusions else \
        {"tic": set(), "gaia": set(), "reason": {}}
    excluded = [i for i in rows if str(st.at[i, "ID"]) in excl["tic"]
                or (pd.notna(gaia_ids[i]) and int(gaia_ids[i]) in excl["gaia"])]
    st.loc[excluded, "tdepth"] = 0.0  # still dilute the aperture, but no scenarios of their own
    fpp_excluded = [str(st.at[i, "ID"]) for i in excluded]
    if fpp_excluded:
        log.info("TIC %d: neighbours %s excluded from the scenarios (pixel fit)", tic_id, fpp_excluded)
    flux = 1 - (1 - g["flux"].to_numpy()) * fr_target
    err = 1.4826 * np.median(np.abs(flux[oot] - np.median(flux[oot])))
    # Without a TRILEGAL population the background scenarios cannot be computed and are dropped.
    # The FPP is then OPTIMISTIC (unresolved background EBs not considered) and is labelled so.
    dropped = []
    if trilegal is None:
        dropped = ["DTP", "DEB", "BTP", "BEB"]
        log.warning("TIC %d: no TRILEGAL population, background scenarios %s dropped (FPP optimistic)",
                    tic_id, dropped)
    runs, all_probs = [], []
    for k in range(n_runs):
        np.random.seed(seed + k)
        target.calc_probs(time=dt[sel], flux_0=flux[sel], flux_err_0=float(err), P_orb=float(P), N=n_draws,
                          exptime=float(g["cadence_s"].iloc[0]) / 86400, verbose=0, drop_scenario=dropped)
        pr = target.probs[["scenario", "ID", "prob"]].copy()
        pr["ID"] = pr["ID"].astype(str)
        all_probs.append(pr.assign(run=k))
        runs.append({"FPP": float(target.FPP), "NFPP": float(target.NFPP), **scenario_groups(pr)})
    runs = pd.DataFrame(runs)
    probs = pd.concat(all_probs).groupby(["scenario", "ID"], as_index=False)["prob"].mean()
    return {"FPP": float(runs["FPP"].mean()), "FPP_std": float(runs["FPP"].std(ddof=1)) if n_runs > 1 else None,
            "NFPP": float(runs["NFPP"].mean()), "NFPP_std": float(runs["NFPP"].std(ddof=1)) if n_runs > 1 else None,
            "FPP_runs": runs["FPP"].round(4).tolist(), "n_runs": n_runs,
            **{f"prob_{c}": float(runs[c].mean()) for c in SCENARIO_GROUPS}, "n_draws": n_draws,
            "dropped_scenarios": dropped, "fpp_background": trilegal is not None,
            "fpp_method": "triceratops + TRILEGAL 1.6" if trilegal else "triceratops, no background (optimistic)",
            "trilegal_error": trilegal_error, "aperture_pixels": int(len(ap)), "target_flux_ratio": fr_target,
            "n_neighbour_scenarios": len(rows) - len(excluded), "neighbour_parameters": nb_log,
            "n_neighbours_filled": sum(any(p in r for p in ("Teff", "rad", "mass")) for r in nb_log),
            "n_neighbours_solar_default": sum("unfilled" in r for r in nb_log), "neighbour_error": nb_error,
            "tic_duplicates_removed": duplicates, "fpp_excluded": fpp_excluded,
            "fpp_exclusion_reason": {k: v for k, v in excl["reason"].items() if k in fpp_excluded},
            "scenarios": probs.sort_values("prob", ascending=False).head(8).to_dict("records")}


# TRICERATOPS scenarios by where the signal comes from
SCENARIO_GROUPS = {
    "target_planet": ("TP", "PTP", "DTP"),                          # a planet on the target
    "target_eb": ("EB", "EBx2P", "PEB", "PEBx2P", "DEB", "DEBx2P"),  # the target is an eclipsing binary
    "companion": ("STP", "SEB", "SEBx2P"),                           # a bound, unresolved companion
    "neighbour": ("NTP", "NEB", "NEBx2P"),                           # a catalogued star in the aperture
    "background": ("BTP", "BEB", "BEBx2P"),                          # an uncatalogued star behind the target
}


def scenario_groups(probs: pd.DataFrame) -> dict:
    """Sum the scenario probabilities by origin (keys of SCENARIO_GROUPS)."""
    return {g: float(probs.loc[probs["scenario"].isin(sc), "prob"].sum()) for g, sc in SCENARIO_GROUPS.items()}


# ------------------------------------------------------------------------------------------------
# driver
# ------------------------------------------------------------------------------------------------


# Bump when the follow-up steps or their stored fields change: summaries written by an older version are
# not "up to date" and are redone by followup_sector() (resume skips only current ones).
FOLLOWUP_VERSION = 3  # 3: difference_image.status from pseudo-transit null images
# which vetting survivors get followed up: the signals not explained by a known object
FOLLOWUP_CLASSES = ("new", "new_in_known_system")


def pixel_verdict(di: dict) -> dict:
    """Read the difference image: on target, off target (conclusive or weak), inconclusive, or missing.

    Difference images with a `status` (on_target / off_target / inconclusive, from difference_image_stats)
    use it directly: its thresholds come from pseudo-transit null images, so an off-target loss there beats
    what noise produces in that cutout, and "off" is strong enough to reject a signal on its own.
    Older summaries without `status` fall back to fixed thresholds:
    "inconclusive": no significant flux loss anywhere (strongest < 5 sigma, < 3 sigma at the target).
    "off": the CENTROID of the loss is > 3 sigma and > 0.5 px from both the light centroid and the catalog
    position. "off_weak": off_target only because the strongest single pixel sits elsewhere. With few transits
    and per-pixel noise from the transit-to-transit scatter this can be a noise pixel (TIC 372868160: S/N 12 at
    10 px, centroid offset 0.3 sigma), so it counts only after the signal has been shown to repeat.
    """
    if not di or "off_target" not in di:
        return {"pixels": "missing" if not di else "error", "offset": None}
    off = f"{di['offset_px']:.2f} ± {di['offset_err_px']:.2f} px"
    status = di.get("status")
    if status is not None:
        return {"pixels": {"on_target": "on", "off_target": "off", "inconclusive": "inconclusive"}[status],
                "offset": off}
    if di["strongest_loss_snr"] < 5 and di["loss_snr_at_target"] < 3:
        return {"pixels": "inconclusive", "offset": off}
    if not di["off_target"]:
        return {"pixels": "on", "offset": off}
    conclusive = di.get("offset_sigma", 0) > 3 and min(di["offset_px"], di.get("offset_from_catalog_px", np.inf)) > 0.5
    return {"pixels": "off" if conclusive else "off_weak", "offset": off}


# verdict categories, best first; the report ranks candidates in this order
CATEGORIES = ["validated", "candidate", "ambiguous", "incomplete", "untestable", "error",
              "not_repeating", "eb", "blend"]
REJECTED = {"not_repeating", "eb", "blend"}


def verdict(summary: dict) -> dict:
    """The follow-up verdict of one candidate from its summary.json: the first check that explains it.

    Order: (1) flux loss conclusively off target -> blend (an EB on a neighbour if it is also deeper at 2P);
    (2) no sector could test repetition -> untestable; (3) deeper at 2P and not repeating at P -> EB found at
    half its period; (4) does not repeat -> not_repeating; (5) weak off-target pixels -> blend; (6) fitted
    stellar density outside 1/3-3x the catalogue -> blend or wrong host; (7) an FPP computed without the
    background scenarios -> incomplete (it is optimistic); (8) FPP < 0.015 and NFPP < 1e-3 -> validated-level,
    FPP < 0.5 -> candidate (CTOI-worthy), else ambiguous; (9) fit or FPP missing -> incomplete.
    The FPP used is the one with the pixel-based neighbour exclusions when there are any (fpp_excluded);
    FPP_without_exclusions is reported next to it.
    Repetition is recomputed from the stored sector table, so rule changes apply without new searches.
    """
    if "error" in summary and "sector_search" not in summary:
        return {"category": "error", "verdict": f"follow-up failed: {summary['error']}"}
    sec = pd.DataFrame(summary["sector_search"])
    rp = repetition(sec, summary["P"], summary["depth_ppm"], summary["discovery_sector"])
    pv = pixel_verdict(summary.get("difference_image", {}))
    out = {"repeats": bool(rp["repeats"]) if rp["testable"] else None,
           "n_same_period": rp["n_same_period"], "n_other_sectors": rp["n_other_sectors"],
           "n_informative": rp["n_informative"], "repetition_basis": rp["basis"],
           "n_deeper_at_2P": rp["n_alias_period"], "pixels": pv["pixels"], "offset": pv["offset"],
           "density_ratio": np.nan, "rp_fit_earth": np.nan, "FPP": np.nan, "NFPP": np.nan,
           "FPP_std": np.nan, "fpp_background": None, "n_fpp_excluded": 0, "FPP_without_exclusions": np.nan}
    f = summary.get("fit") or {}
    if "rho_ratio_median" in f:
        out["density_ratio"] = float(f["rho_ratio_median"])
        out["rp_fit_earth"] = float(f["rp_earth"]["median"])
    fp = summary.get("fpp") or {}
    if "FPP" in fp:
        out["FPP"], out["NFPP"] = float(fp["FPP"]), float(fp["NFPP"])  # unrounded: 0.0148 must pass "< 0.015"
        out["FPP_std"] = float(fp.get("FPP_std", np.nan))
        out["fpp_background"] = bool(fp.get("fpp_background", False))
        out["n_fpp_excluded"] = len(fp.get("fpp_excluded") or [])
        out["FPP_without_exclusions"] = float((fp.get("without_exclusions") or {}).get("FPP", np.nan))
    density_ok = None if np.isnan(out["density_ratio"]) else bool(1 / 3 < out["density_ratio"] < 3)

    if pv["pixels"] == "off":
        cat, text = ("eb", "eclipsing binary on a neighbour, found at half its period") if rp["repeats_at_alias"] \
            else ("blend", "flux loss off target: blend on a neighbour")
    elif not rp["testable"]:
        cat, text = "untestable", "no verdict: no other sector could show it"
    elif not rp["repeats"] and rp["repeats_at_alias"]:
        cat, text = "eb", "eclipsing binary at half its period (deeper at 2P)"
    elif not rp["repeats"]:
        cat, text = "not_repeating", "not confirmed: does not repeat"
    elif pv["pixels"] == "off_weak":
        cat, text = "blend", "strongest flux loss off target: blend"
    elif density_ok is False:
        cat, text = "blend", "inconsistent stellar density: blend or wrong host"
    elif np.isfinite(out["FPP"]) and not out["fpp_background"]:
        # without the background scenarios (no TRILEGAL population) the FPP is optimistic: never a candidate
        cat, text = "incomplete", "FPP without background scenarios (optimistic)"
    elif np.isfinite(out["FPP"]) and out["FPP"] < 0.015 and out["NFPP"] < 1e-3:
        cat, text = "validated", "validated-level FPP"
    elif np.isfinite(out["FPP"]) and out["FPP"] < 0.5:
        cat, text = "candidate", "planet candidate (CTOI-worthy)"
    elif np.isfinite(out["FPP"]):
        cat, text = "ambiguous", "ambiguous (FPP >= 0.5)"
    else:
        cat, text = "incomplete", "incomplete: fit or FPP missing" + (
            f" ({fp['error']})" if "error" in fp else "")
    return {"category": cat, "verdict": text, **out}


def _cheap_reject(summary: dict) -> bool:
    """True when repetition and pixels already explain the signal, so the fit and FPP are not needed."""
    v = verdict({k: summary[k] for k in summary if k not in ("fit", "fpp")})
    return v["category"] in REJECTED or v["category"] == "untestable"


def _reusable(old: dict, key: str, P: float, T0: float) -> bool:
    """An earlier fit / FPP result for the same ephemeris that can be kept instead of recomputed (both are
    minutes of MCMC / TRICERATOPS and do not change unless the signal does). An FPP computed without the
    background scenarios (no `fpp_background`, the pre-TRILEGAL runs) is not reused."""
    r = old.get(key) or {}
    if not r or "error" in r or abs(old.get("P", np.nan) - P) > 1e-6 or abs(old.get("T0", np.nan) - T0) > 1e-6:
        return False
    return key != "fpp" or bool(r.get("fpp_background"))


def followup_one(tic_id: int, sector: int, skip_fit: bool = False, skip_fpp: bool = False,
                 threads: int = 1, auto: bool = False, reuse: bool = True) -> dict:
    """All follow-up checks for one candidate; writes data/followup/{tic}/summary.json.

    auto=True runs the expensive fit and FPP only when repetition and pixels have not already rejected the
    signal (the pipeline's mode). skip_fit / skip_fpp skip them regardless. reuse=True keeps an earlier fit /
    FPP for the same ephemeris (see _reusable); the repetition searches and the difference image are always
    redone, since they are where the follow-up rules change.
    """
    known = pd.read_parquet(DATA_DIR / "search" / f"sector_{sector}_known.parquet").set_index("tic_id")
    targets = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}.parquet").set_index("tic_id")
    s, star = known.loc[tic_id], targets.loc[tic_id]
    out_dir = FOLLOWUP_DIR / str(tic_id)
    out_dir.mkdir(parents=True, exist_ok=True)
    # keep the results of steps skipped in this run (--skip-fit / --skip-fpp re-run only the cheap checks)
    old_path = out_dir / "summary.json"
    summary = json.loads(old_path.read_text()) if old_path.exists() else {}
    summary.pop("error", None)
    if reuse:
        skip_fit = skip_fit or _reusable(summary, "fit", float(s.period), float(s.T0))
        skip_fpp = skip_fpp or _reusable(summary, "fpp", float(s.period), float(s.T0))
    summary.update({"tic_id": tic_id, "discovery_sector": sector, "P": float(s.period), "T0": float(s.T0),
                    "depth_ppm": float(s.depth_ppm), "duration_h": float(s.duration_h),
                    "followup_version": FOLLOWUP_VERSION})

    ms = fetch_sectors(tic_id, out_dir)
    summary["sectors"] = sorted(int(x) for x in ms["sector"].unique())
    sec = sector_search(ms, s.period, star.rad, star.mass, sector, threads=threads, dur_h=s.duration_h,
                        depth_ppm=s.depth_ppm)
    sec = pd.concat([sec, alias_search(ms, s.period, star.rad, star.mass, threads=threads)], ignore_index=True)
    sec.to_csv(out_dir / "sectors.csv", index=False)
    summary["sector_search"] = sec.to_dict("records")
    summary["repetition"] = repetition(sec, s.period, s.depth_ppm, sector)
    comb = sec[sec["data"] == "all sectors"].iloc[0]
    P_best, T0_best = float(comb["period"]), float(comb["T0"])

    try:
        di = difference_image(tic_id, star, s.period, s.T0, s.duration_h, sector, out_dir)
        di.pop("_tpf_wcs", None)
        summary["difference_image"] = di
    except Exception as ex:
        summary["difference_image"] = {"error": f"{type(ex).__name__}: {ex}"}
    try:
        g = gaia_neighbours(target_coord(star, s.T0), s.depth_ppm)
        g.to_csv(out_dir / "gaia.csv", index=False)
        summary["gaia_could_cause_dip"] = int(g["could_cause_dip"].iloc[1:].sum()) if len(g) else 0
    except Exception as ex:
        summary["gaia_could_cause_dip"] = None
        log.warning("TIC %d: Gaia query failed: %s", tic_id, ex)

    if auto and _cheap_reject(summary):
        skip_fit = skip_fpp = True
        summary["full_chain"] = False
        log.info("TIC %d: rejected by repetition/pixels, fit and FPP skipped", tic_id)
    else:
        summary["full_chain"] = not (skip_fit and skip_fpp) or ("fit" in summary and "fpp" in summary)
    if not skip_fit:
        try:
            fit, chain = fit_transit(ms, P_best, T0_best, s.duration_h, s.depth_ppm, star.rad, star.mass)
            np.savez_compressed(out_dir / "fit_samples.npz", chain=chain)
            (out_dir / "fit.json").write_text(json.dumps(fit, indent=2))
            summary["fit"] = fit
        except Exception as ex:
            summary["fit"] = {"error": f"{type(ex).__name__}: {ex}"}
    if not skip_fpp:
        try:
            summary["fpp"] = fpp(tic_id, ms, s.period, s.T0, s.duration_h, s.depth_ppm, sector)
        except Exception as ex:
            summary["fpp"] = {"error": f"{type(ex).__name__}: {ex}"}
    # a step skipped here may have been (re)computed by another process meanwhile: keep the newer result on disk
    if old_path.exists():
        now = json.loads(old_path.read_text())
        for key, skipped in (("fit", skip_fit), ("fpp", skip_fpp)):
            if skipped and key in now:
                summary[key] = now[key]
    summary["verdict"] = verdict(summary)
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    log.info("TIC %d: %s (%s)", tic_id, summary["verdict"]["verdict"], out_dir / "summary.json")
    return summary


def survivors(known: pd.DataFrame) -> pd.DataFrame:
    """Vetting survivors that no known object explains: the signals the follow-up stage checks."""
    return known[known["passes_all"] & known["known_class"].isin(FOLLOWUP_CLASSES)]


def up_to_date(tic_id: int, row: pd.Series) -> bool:
    """A summary.json exists, was written by this FOLLOWUP_VERSION, for the same ephemeris, and has a verdict."""
    path = FOLLOWUP_DIR / str(tic_id) / "summary.json"
    if not path.exists():
        return False
    s = json.loads(path.read_text())
    return (s.get("followup_version") == FOLLOWUP_VERSION and "verdict" in s and "error" not in s
            and abs(s.get("P", np.nan) - row["period"]) < 1e-6 and abs(s.get("T0", np.nan) - row["T0"]) < 1e-6)


def _is_network_error(ex: Exception) -> bool:
    """MAST / TESScut / Vizier hiccups (worth a retry), as opposed to a bug or bad data."""
    names = {c.__name__ for c in type(ex).__mro__}
    return bool(names & {"RemoteServiceError", "ConnectionError", "Timeout", "HTTPError", "ReadTimeout",
                         "TimeoutError", "URLError", "SearchError"}) or "Request failed" in str(ex)


def _safe_one(tic_id: int, sector: int, auto: bool, retries: int = 3, wait_s: float = 30.0) -> dict:
    """followup_one that never raises: network errors are retried, a failure is recorded in summary.json and
    the run goes on."""
    import warnings as _w

    logging.getLogger("lightkurve").setLevel(logging.ERROR)
    _w.filterwarnings("ignore")
    try:
        for attempt in range(1, retries + 1):
            try:
                return followup_one(tic_id, sector, auto=auto)
            except Exception as ex:
                if attempt == retries or not _is_network_error(ex):
                    raise
                log.warning("TIC %d: %s, retry %d/%d in %.0f s", tic_id, ex, attempt, retries - 1, wait_s * attempt)
                time.sleep(wait_s * attempt)
    except Exception as ex:
        log.exception("TIC %d: follow-up failed", tic_id)
        out_dir = FOLLOWUP_DIR / str(tic_id)
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "summary.json"
        s = json.loads(path.read_text()) if path.exists() else {}
        s.update({"tic_id": tic_id, "discovery_sector": sector, "followup_version": FOLLOWUP_VERSION,
                  "error": f"{type(ex).__name__}: {ex}"})
        s.pop("verdict", None)
        path.write_text(json.dumps(s, indent=2, default=float))
        return s


def followup_table(summaries: list[dict]) -> pd.DataFrame:
    """One row per candidate: verdict and the key numbers behind it."""
    rows = []
    for s in summaries:
        v = verdict(s) if "error" not in s else {"category": "error", "verdict": f"follow-up failed: {s['error']}"}
        rows.append({"tic_id": int(s["tic_id"]), "n_sectors": len(s.get("sectors", [])), **v,
                     "gaia_could_cause_dip": s.get("gaia_could_cause_dip"), "full_chain": s.get("full_chain"),
                     "error": s.get("error")})
    out = pd.DataFrame(rows)
    if len(out):
        out["rank_key"] = out["category"].map(CATEGORIES.index)
        out = out.sort_values(["rank_key", "tic_id"]).drop(columns="rank_key").reset_index(drop=True)
    return out


def refresh_verdicts(sector: int) -> pd.DataFrame:
    """Recompute every survivor's verdict from its stored summary.json (no searches, no network) and rewrite
    the sector table: for when the verdict rules change or another process has updated a fit or FPP."""
    known = pd.read_parquet(DATA_DIR / "search" / f"sector_{sector}_known.parquet")
    summaries = []
    for t in survivors(known)["tic_id"]:
        path = FOLLOWUP_DIR / str(int(t)) / "summary.json"
        if not path.exists():
            summaries.append({"tic_id": int(t), "error": "not followed up yet"})
            continue
        s = json.loads(path.read_text())
        if "error" not in s:
            s["verdict"] = verdict(s)
            path.write_text(json.dumps(s, indent=2, default=float))
        summaries.append(s)
    table = followup_table(summaries)
    table.to_parquet(FOLLOWUP_DIR / f"sector_{sector}_followup.parquet", index=False)
    return table


def default_n_jobs() -> int:
    """At most half the cores: each worker runs TLS, and the machine is shared with other searches."""
    return max(1, (os.cpu_count() or 2) // 2)


def followup_sector(sector: int, n_jobs: int | None = None, force: bool = False, tics: list[int] | None = None,
                    auto: bool = True) -> pd.DataFrame:
    """Follow up every survivor of `sector` (or `tics`), in parallel and resumably.

    Writes data/followup/sector_{N}_followup.parquet. Candidates whose summary is up to date are not redone
    (force=True redoes them); a candidate that fails is recorded with its error and does not stop the run.
    """
    import joblib

    known = pd.read_parquet(DATA_DIR / "search" / f"sector_{sector}_known.parquet").set_index("tic_id")
    todo_tab = known.loc[tics] if tics else survivors(known.reset_index()).set_index("tic_id")
    redo = [int(t) for t, r in todo_tab.iterrows() if force or not up_to_date(int(t), r)]
    log.info("sector %d: %d candidates to follow up, %d up to date, %d to run", sector, len(todo_tab),
             len(todo_tab) - len(redo), len(redo))
    if redo:
        n_jobs = min(n_jobs or default_n_jobs(), len(redo))
        joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(_safe_one)(t, sector, auto) for t in redo)
    summaries = [json.loads((FOLLOWUP_DIR / str(int(t)) / "summary.json").read_text()) for t in todo_tab.index]
    table = followup_table(summaries)
    FOLLOWUP_DIR.mkdir(parents=True, exist_ok=True)
    if not tics:
        path = FOLLOWUP_DIR / f"sector_{sector}_followup.parquet"
        table.to_parquet(path, index=False)
        log.info("wrote %s: %s", path, table["category"].value_counts().to_dict() if len(table) else {})
    return table


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sector", type=int, default=int(os.environ.get("SECTOR", 0)) or None)
    p.add_argument("--tics", default="", help="comma-separated TIC IDs (default: every vetting survivor of the "
                                              "sector that no known object explains, see survivors())")
    p.add_argument("--skip-fit", action="store_true")
    p.add_argument("--skip-fpp", action="store_true")
    p.add_argument("--full-chain", action="store_true",
                   help="sector mode: fit and FPP even for signals already rejected by repetition/pixels")
    p.add_argument("--n-jobs", type=int, default=None, help="parallel candidates (default: half the cores)")
    p.add_argument("--force", action="store_true", help="redo candidates whose summary is up to date")
    p.add_argument("--threads", type=int, default=1, help="TLS threads for the multi-sector searches")
    p.add_argument("--verdict-only", action="store_true",
                   help="only recompute the verdicts and the sector table from the stored summaries")
    a = p.parse_args()
    if a.sector is None:
        p.error("no --sector given and SECTOR not set in .env")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.getLogger("lightkurve").setLevel(logging.ERROR)
    tics = [int(x) for x in a.tics.split(",") if x.strip()]
    if tics and (a.skip_fit or a.skip_fpp):  # hand-run quick checks on chosen stars, as before
        for tic in tics:
            followup_one(tic, a.sector, a.skip_fit, a.skip_fpp, a.threads)
        return
    if a.verdict_only:
        table = refresh_verdicts(a.sector)
    else:
        table = followup_sector(a.sector, a.n_jobs, a.force, tics or None, auto=not a.full_chain)
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(table[["tic_id", "category", "verdict", "repeats", "n_same_period", "n_informative", "pixels",
                     "density_ratio", "FPP"]].to_string() if len(table) else "no candidates")


if __name__ == "__main__":
    main()
