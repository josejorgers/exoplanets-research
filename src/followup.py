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


def difference_image(tic_id: int, star: pd.Series, P: float, T0: float, dur_h: float, sector: int,
                     out_dir: Path, size: int = 15, n_boot: int = 200, seed: int = 0) -> dict:
    """Mean out-of-transit minus mean in-transit image, transit by transit, from a TESScut cutout.

    The difference image shows where the flux disappears during transit. Its centroid is compared
    with the target's (proper-motion corrected) position; the uncertainty comes from bootstrapping
    over transits. An offset of more than ~1 pixel at high significance means the source is a
    neighbour.
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
    dur = dur_h / 24
    epoch = np.round((t - T0) / P).astype(int)
    dt = t - (T0 + epoch * P)
    diffs, oots = [], []
    for k in np.unique(epoch):
        i_in = (epoch == k) & (np.abs(dt) < 0.5 * dur)  # the whole transit: all the signal there is
        i_oot = (epoch == k) & (np.abs(dt) > 0.75 * dur) & (np.abs(dt) < 2.5 * dur)
        if i_in.sum() >= 3 and i_oot.sum() >= 6:
            oot = np.median(cube[i_oot], axis=0)
            diffs.append(oot - np.mean(cube[i_in], axis=0))
            oots.append(oot)
    if len(diffs) < 2:
        return {"n_transits": len(diffs), "error": "fewer than 2 usable transits in the cutout"}
    diffs, oots = np.array(diffs), np.array(oots)
    diff, oot = diffs.mean(axis=0), oots.mean(axis=0)
    tgt = target_coord(star, float(np.median(t)))
    tx, ty = tpf.wcs.world_to_pixel(tgt)

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

    cx, cy = centroid(diff)
    ox, oy = centroid(oot)
    rng = np.random.default_rng(seed)
    boot = np.array([centroid(diffs[rng.integers(0, len(diffs), len(diffs))].mean(axis=0)) for _ in range(n_boot)])
    boot = boot[np.all(np.isfinite(boot), axis=1)]
    sx, sy = boot.std(axis=0) if len(boot) > 10 else (np.nan, np.nan)
    sig = float(np.hypot(sx, sy))
    off = float(np.hypot(cx - ox, cy - oy))  # difference centroid vs the target's own light centroid
    off_cat = float(np.hypot(cx - tx, cy - ty))  # vs the catalog position
    # where is the strongest flux loss anywhere in the cutout, in units of its scatter across transits?
    snr_img = diff / (diffs.std(axis=0) / np.sqrt(len(diffs)))
    gy, gx = np.unravel_index(np.nanargmax(snr_img), snr_img.shape)
    ix, iy = int(round(float(tx))), int(round(float(ty)))
    res = {
        "sector": sector, "n_transits": len(diffs), "target_px": [float(tx), float(ty)],
        "diff_centroid_px": [cx, cy], "oot_centroid_px": [ox, oy],
        "offset_px": off, "offset_arcsec": off * PIXEL_ARCSEC, "offset_err_px": sig,
        "offset_sigma": off / sig if sig > 0 else np.nan, "offset_from_catalog_px": off_cat,
        "strongest_loss_px": [int(gx), int(gy)], "strongest_loss_snr": float(snr_img[gy, gx]),
        "strongest_loss_dist_px": float(np.hypot(gx - tx, gy - ty)),
        "loss_snr_at_target": float(snr_img[iy, ix]) if 0 <= ix < size and 0 <= iy < size else np.nan,
    }
    sky = tpf.wcs.pixel_to_world(gx, gy)
    res["strongest_loss_radec"] = [float(sky.ra.deg), float(sky.dec.deg)]
    # off target: the centroid of the flux loss is significantly away from the star's own light AND from its
    # catalog position, or the strongest loss in the cutout sits > 1.5 px away and is clearly stronger than at
    # the target. Both positions are needed: in a crowded field the out-of-transit light centroid is pulled
    # towards the neighbours (TIC 214311908: 0.76 px from the light centroid, 0.30 px from the catalog position).
    res["off_target"] = bool((res["offset_sigma"] > 3 and min(off, off_cat) > 0.5)
                             or (res["strongest_loss_dist_px"] > 1.5 and res["strongest_loss_snr"] > 5
                                 and res["loss_snr_at_target"] < 0.5 * res["strongest_loss_snr"]))
    np.savez_compressed(out_dir / f"diffimage_s{sector}.npz", diff=diff, oot=oot, target=[tx, ty],
                        centroid=[cx, cy], oot_centroid=[ox, oy], strongest=[gx, gy], snr_img=snr_img, boot=boot)
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
        n_draws: int = 20_000, seed: int = 0) -> dict:
    """TRICERATOPS false-positive probability from the discovery sector's folded light curve."""
    from src import _compat  # noqa: F401  (scipy/pkg_resources shims needed by pytransit)
    import triceratops.triceratops as tr

    np.random.seed(seed)
    g = ms[ms["sector"] == sector]
    dt = (g["time"].to_numpy() - T0 + P / 2) % P - P / 2
    sel = np.abs(dt) < 3 * dur_h / 24
    oot = np.abs(dt) > dur_h / 24
    target = tr.target(ID=tic_id, sectors=np.array([sector]))
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
    target.calc_depths(tdepth=0.5 * depth_ppm / 1e6, all_ap_pixels=[ap])
    fr_target = float(target.stars["fluxratio"].iloc[0])
    target.calc_depths(tdepth=depth_ppm / 1e6 * fr_target, all_ap_pixels=[ap])
    if not (target.stars["tdepth"] > 0).any():
        raise RuntimeError("no star can produce the transit depth; check the depth units")
    flux = 1 - (1 - g["flux"].to_numpy()) * fr_target
    err = 1.4826 * np.median(np.abs(flux[oot] - np.median(flux[oot])))
    # Background-star scenarios need a TRILEGAL Galaxy-model population. When the TRILEGAL server is
    # unreachable this triceratops version announces it will skip them but then crashes in them, so
    # drop them explicitly. The FPP is then OPTIMISTIC (unresolved background EBs not considered).
    dropped = []
    if getattr(target, "trilegal_url", None) is None and getattr(target, "trilegal_fname", None) is None:
        dropped = ["DTP", "DEB", "BTP", "BEB"]
        log.warning("TIC %d: TRILEGAL unavailable, background scenarios %s dropped (FPP optimistic)", tic_id, dropped)
    target.calc_probs(time=dt[sel], flux_0=flux[sel], flux_err_0=float(err), P_orb=float(P), N=n_draws,
                      exptime=float(g["cadence_s"].iloc[0]) / 86400, verbose=0, drop_scenario=dropped)
    probs = target.probs[["scenario", "ID", "prob"]].copy()
    probs["ID"] = probs["ID"].astype(str)
    return {"FPP": float(target.FPP), "NFPP": float(target.NFPP), "n_draws": n_draws,
            "dropped_scenarios": dropped, "aperture_pixels": int(len(ap)), "target_flux_ratio": fr_target,
            "scenarios": probs.sort_values("prob", ascending=False).head(8).to_dict("records")}


# ------------------------------------------------------------------------------------------------
# driver
# ------------------------------------------------------------------------------------------------


def followup_one(tic_id: int, sector: int, skip_fit: bool = False, skip_fpp: bool = False,
                 threads: int = 1) -> dict:
    known = pd.read_parquet(DATA_DIR / "search" / f"sector_{sector}_known.parquet").set_index("tic_id")
    targets = pd.read_parquet(DATA_DIR / "targets" / f"sector_{sector}.parquet").set_index("tic_id")
    s, star = known.loc[tic_id], targets.loc[tic_id]
    out_dir = FOLLOWUP_DIR / str(tic_id)
    out_dir.mkdir(parents=True, exist_ok=True)
    # keep the results of steps skipped in this run (--skip-fit / --skip-fpp re-run only the cheap checks)
    old_path = out_dir / "summary.json"
    summary = json.loads(old_path.read_text()) if old_path.exists() else {}
    summary.update({"tic_id": tic_id, "discovery_sector": sector, "P": float(s.period), "T0": float(s.T0),
                    "depth_ppm": float(s.depth_ppm), "duration_h": float(s.duration_h)})

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

    if not skip_fit:
        fit, chain = fit_transit(ms, P_best, T0_best, s.duration_h, s.depth_ppm, star.rad, star.mass)
        np.savez_compressed(out_dir / "fit_samples.npz", chain=chain)
        (out_dir / "fit.json").write_text(json.dumps(fit, indent=2))
        summary["fit"] = fit
    try:
        di = difference_image(tic_id, star, s.period, s.T0, s.duration_h, sector, out_dir)
        di.pop("_tpf_wcs", None)
        summary["difference_image"] = di
    except Exception as ex:
        summary["difference_image"] = {"error": f"{type(ex).__name__}: {ex}"}
    g = gaia_neighbours(target_coord(star, s.T0), s.depth_ppm)
    g.to_csv(out_dir / "gaia.csv", index=False)
    summary["gaia_could_cause_dip"] = int(g["could_cause_dip"].iloc[1:].sum()) if len(g) else 0
    if not skip_fpp:
        try:
            summary["fpp"] = fpp(tic_id, ms, s.period, s.T0, s.duration_h, s.depth_ppm, sector)
        except Exception as ex:
            summary["fpp"] = {"error": f"{type(ex).__name__}: {ex}"}
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    log.info("TIC %d: wrote %s", tic_id, out_dir / "summary.json")
    return summary


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sector", type=int, default=int(os.environ.get("SECTOR", 0)) or None)
    p.add_argument("--tics", required=True, help="comma-separated TIC IDs")
    p.add_argument("--skip-fit", action="store_true")
    p.add_argument("--skip-fpp", action="store_true")
    p.add_argument("--threads", type=int, default=1, help="TLS threads for the multi-sector searches")
    a = p.parse_args()
    if a.sector is None:
        p.error("no --sector given and SECTOR not set in .env")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.getLogger("lightkurve").setLevel(logging.ERROR)
    for tic in [int(x) for x in a.tics.split(",") if x.strip()]:
        followup_one(tic, a.sector, a.skip_fit, a.skip_fpp, a.threads)


if __name__ == "__main__":
    main()
