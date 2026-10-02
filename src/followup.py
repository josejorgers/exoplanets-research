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


def _tls_window(lc: pd.DataFrame, P: float, rad: float, mass: float, window: float, threads: int = 1):
    t, y, dy = cleaned_array(lc["time"].values, lc["flux"].values, lc["flux_err"].values)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return transitleastsquares(t, y, dy).power(
            R_star=rad, M_star=mass, R_star_min=0.1, R_star_max=1.0, M_star_min=0.05, M_star_max=1.0,
            period_min=P * (1 - window), period_max=P * (1 + window), oversampling_factor=5,
            duration_grid_step=1.05, use_threads=threads, show_progress_bar=False, verbose=False,
        )


def sector_search(ms: pd.DataFrame, P: float, rad: float, mass: float, discovery_sector: int,
                  window: float = 0.02, threads: int = 1) -> pd.DataFrame:
    """TLS restricted to P*(1 +- window): per sector, all but the discovery sector, and all sectors.

    In a narrow window the look-elsewhere effect is small, so a real transit shows up with a high
    SDE in every sector where it has enough transits; noise or a one-sector systematic does not.
    """
    groups = [(f"s{s:03d}", g) for s, g in ms.groupby("sector")]
    groups += [("all but discovery", ms[ms["sector"] != discovery_sector]), ("all sectors", ms)]
    rows = []
    for name, g in groups:
        if len(g) < 500:
            continue
        r = _tls_window(g, P, rad, mass, window, threads)
        rows.append({
            "data": name, "n_points": len(g), "author": "/".join(sorted(g["author"].unique())),
            "SDE": float(r.SDE), "period": float(r.period), "period_unc": float(r.period_uncertainty),
            "T0": float(r.T0), "depth_ppm": (1 - float(r.depth)) * 1e6, "snr": float(r.snr),
            "duration_h": float(r.duration) * 24, "n_transits": int(r.distinct_transit_count),
            "odd_even_mismatch": float(r.odd_even_mismatch),
        })
        log.info("%s: SDE %.1f P %.5f depth %.0f ppm (%d transits)", name, r.SDE, r.period, (1 - r.depth) * 1e6,
                 r.distinct_transit_count)
    return pd.DataFrame(rows)


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
                nwalkers: int = 32, nsteps: int = 3000, seed: int = 0) -> tuple[dict, np.ndarray]:
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
    dur = dur_h / 24
    epoch = np.round((t - T0) / P).astype(int)
    dt = t - (T0 + epoch * P)
    diffs, oots = [], []
    for k in np.unique(epoch):
        i_in = (epoch == k) & (np.abs(dt) < 0.4 * dur)
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

    def centroid(img):
        """Flux-weighted centroid of the positive pixels within 2 px of the image's brightest pixel."""
        yy, xx = np.indices(img.shape)
        py, px = np.unravel_index(np.nanargmax(img), img.shape)
        w = np.where((np.hypot(xx - px, yy - py) <= 2) & (img > 0), img, 0.0)
        return float((w * xx).sum() / w.sum()), float((w * yy).sum() / w.sum())

    cx, cy = centroid(diff)
    ox, oy = centroid(oot)
    rng = np.random.default_rng(seed)
    boot = np.array([centroid(diffs[rng.integers(0, len(diffs), len(diffs))].mean(axis=0)) for _ in range(n_boot)])
    sx, sy = boot.std(axis=0)
    off = float(np.hypot(cx - tx, cy - ty))
    sig = float(np.hypot(sx, sy))
    # significance of the difference image at the target: depth / its scatter across transits
    ix, iy = int(round(float(tx))), int(round(float(ty)))
    res = {
        "sector": sector, "n_transits": len(diffs), "target_px": [float(tx), float(ty)],
        "diff_centroid_px": [cx, cy], "oot_centroid_px": [ox, oy],
        "offset_px": off, "offset_arcsec": off * PIXEL_ARCSEC, "offset_err_px": sig,
        "offset_sigma": off / sig if sig > 0 else np.nan,
        "diff_oot_offset_px": float(np.hypot(cx - ox, cy - oy)),
        "depth_at_target_pixel": float(diff[iy, ix] / oot[iy, ix]) if 0 <= ix < size and 0 <= iy < size else np.nan,
    }
    np.savez_compressed(out_dir / f"diffimage_s{sector}.npz", diff=diff, oot=oot, target=[tx, ty],
                        centroid=[cx, cy], boot=boot)
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
        n_draws: int = 100_000, seed: int = 0) -> dict:
    """TRICERATOPS false-positive probability from the discovery sector's folded light curve."""
    from src import _compat  # noqa: F401  (scipy/pkg_resources shims needed by pytransit)
    import triceratops.triceratops as tr

    np.random.seed(seed)
    g = ms[ms["sector"] == sector]
    dt = (g["time"].to_numpy() - T0 + P / 2) % P - P / 2
    sel = np.abs(dt) < 3 * dur_h / 24
    oot = np.abs(dt) > dur_h / 24
    flux = g["flux"].to_numpy()
    err = 1.4826 * np.median(np.abs(flux[oot] - np.median(flux[oot])))
    target = tr.target(ID=tic_id, sectors=np.array([sector]))
    target.calc_depths(tdepth=depth_ppm)
    target.calc_probs(time=dt[sel], flux_0=flux[sel], flux_err_0=float(err), P_orb=float(P), N=n_draws,
                      exptime=float(g["cadence_s"].iloc[0]) / 86400, verbose=0)
    probs = target.probs[["scenario", "ID", "prob"]].copy()
    probs["ID"] = probs["ID"].astype(str)
    return {"FPP": float(target.FPP), "NFPP": float(target.NFPP), "n_draws": n_draws,
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
    summary = {"tic_id": tic_id, "discovery_sector": sector, "P": float(s.period), "T0": float(s.T0),
               "depth_ppm": float(s.depth_ppm), "duration_h": float(s.duration_h)}

    ms = fetch_sectors(tic_id, out_dir)
    summary["sectors"] = sorted(int(x) for x in ms["sector"].unique())
    sec = sector_search(ms, s.period, star.rad, star.mass, sector, threads=threads)
    sec.to_csv(out_dir / "sectors.csv", index=False)
    summary["sector_search"] = sec.to_dict("records")
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
