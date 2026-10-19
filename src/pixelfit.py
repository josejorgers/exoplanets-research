"""Which star in a crowded TESS cutout loses light during the transit? A PRF-based pixel test.

The difference image (out-of-transit minus in-transit, followup.difference_image) shows where flux
disappears. Its centroid alone cannot tell apart stars closer than about a pixel, and in a crowded
field the out-of-transit light centroid is itself pulled by the neighbours. This module does the
full version of the test:

1. Scene model: the out-of-transit image is fitted as the TESS pixel response function (PRF, the
   SPOC model from MAST) placed at every catalogue star, with a global pointing offset. This gives
   each star's flux in the cutout and checks that PRF and positions are right.
2. Source test: the mean difference image is fitted with the PRF of one star at a time (plus a
   constant). The star whose PRF fits best is the source; the chi-square difference between stars
   measures how strongly the others are excluded. Per-pixel errors and the null distribution come
   from "fake transits" at out-of-transit phases of the same cutout.
3. The fitted amplitude divided by the star's scene flux is the eclipse depth that star would need.
   A depth above 100 % rules the star out regardless of the fit.
"""

from __future__ import annotations

import os
import re
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.ndimage import map_coordinates
from scipy.optimize import minimize, nnls

PRF_URL = "https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/{start}/cam{cam}_ccd{ccd}/"
PRF_SUB = 9  # PRF samples per TESS pixel
PRF_GRID = (1, 513, 1025, 1536, 2048)  # CCD rows of the PRF grid points; columns are these + 44


# ------------------------------------------------------------------------------------------------
# PRF
# ------------------------------------------------------------------------------------------------


def _prf_files(start: str, cam: int, ccd: int) -> list[str]:
    html = urllib.request.urlopen(PRF_URL.format(start=start, cam=cam, ccd=ccd), timeout=60).read().decode()
    return sorted(set(re.findall(r'href="([^"]+prf-[^"]+\.fits)"', html)))


def load_prf(cam: int, ccd: int, col: float, row: float, sector: int, cache_dir: Path) -> np.ndarray:
    """Supersampled PRF image (117 x 117: 13 x 13 pixels, 9 samples each) at a CCD position.

    Bilinear interpolation between the four surrounding PRF grid points, as the SPOC README
    prescribes (accurate to 1-5 %). `col`, `row` are tess-point pixel coordinates.
    """
    start = "start_s0001" if sector <= 3 else "start_s0004"
    cache_dir = Path(cache_dir) / start / f"cam{cam}_ccd{ccd}"
    cache_dir.mkdir(parents=True, exist_ok=True)
    names = None
    rows = np.array(PRF_GRID, float)
    cols = rows + 44

    def grid_image(r: int, c: int) -> np.ndarray:
        nonlocal names
        pat = f"row{r:04d}-col{c:04d}.fits"
        hit = list(cache_dir.glob(f"*{pat}"))
        if not hit:
            names = names or _prf_files(start, cam, ccd)
            name = next(n for n in names if n.endswith(pat))
            urllib.request.urlretrieve(PRF_URL.format(start=start, cam=cam, ccd=ccd) + name, cache_dir / name)
            hit = [cache_dir / name]
        return fits.getdata(hit[0], 0).astype(float)

    i = int(np.clip(np.searchsorted(rows, row) - 1, 0, len(rows) - 2))
    j = int(np.clip(np.searchsorted(cols, col) - 1, 0, len(cols) - 2))
    fr = np.clip((row - rows[i]) / (rows[i + 1] - rows[i]), 0, 1)
    fc = np.clip((col - cols[j]) / (cols[j + 1] - cols[j]), 0, 1)
    img = ((1 - fr) * (1 - fc) * grid_image(int(rows[i]), int(cols[j]))
           + (1 - fr) * fc * grid_image(int(rows[i]), int(cols[j + 1]))
           + fr * (1 - fc) * grid_image(int(rows[i + 1]), int(cols[j]))
           + fr * fc * grid_image(int(rows[i + 1]), int(cols[j + 1])))
    return img


def prf_stamp(prf: np.ndarray, x0: float, y0: float, shape: tuple[int, int]) -> np.ndarray:
    """Fraction of a star's flux falling in each pixel of a cutout, for a star at pixel (x0, y0).

    x is the column (axis 1) and y the row (axis 0); pixel centres are at integer positions. The
    supersampled PRF is sampled at the offset of each pixel centre from the star. Pixels farther
    than 6 pixels from the star get nothing (the PRF model spans 13 x 13 pixels).
    """
    c = (prf.shape[0] - 1) / 2
    yy, xx = np.indices(shape, dtype=float)
    u = c + PRF_SUB * (yy - y0)
    v = c + PRF_SUB * (xx - x0)
    out = map_coordinates(prf, [u.ravel(), v.ravel()], order=1, mode="constant", cval=0.0).reshape(shape)
    # normalise with the full 13 x 13 response at the same sub-pixel phase, so the stamp sums to the
    # fraction of the star's light inside the cutout (< 1 for a star near the edge)
    fy, fx = (y0 - np.floor(y0)), (x0 - np.floor(x0))
    k = np.arange(-6, 7, dtype=float)
    uu, vv = np.meshgrid(c + PRF_SUB * (k - fy), c + PRF_SUB * (k - fx), indexing="ij")
    tot = map_coordinates(prf, [uu.ravel(), vv.ravel()], order=1, mode="constant", cval=0.0).sum()
    return out / tot if tot > 0 else out


# ------------------------------------------------------------------------------------------------
# scene model of the out-of-transit image
# ------------------------------------------------------------------------------------------------


def scene_fit(oot: np.ndarray, xy: np.ndarray, tmag: np.ndarray, prf: np.ndarray, n_free: int = 8,
              mask: np.ndarray | None = None) -> dict:
    """Fit the out-of-transit image as PRF x star fluxes + background, with one global pixel offset.

    The `n_free` brightest stars get their own flux; the fainter ones share a single scale on
    10^(-0.4 Tmag). The offset (dx, dy) absorbs WCS and proper-motion errors. Returns the fitted
    flux of every star, the offset, and the model image.
    """
    mask = np.ones(oot.shape, bool) if mask is None else mask
    order = np.argsort(tmag)
    free = order[:n_free]
    rest = order[n_free:]
    w0 = 10 ** (-0.4 * (tmag - tmag.min()))

    def design(dx, dy):
        stamps = np.array([prf_stamp(prf, x + dx, y + dy, oot.shape) for x, y in xy])
        cols = [stamps[k] for k in free]
        cols.append(np.tensordot(w0[rest], stamps[rest], axes=1) if len(rest) else np.zeros(oot.shape))
        cols.append(np.ones(oot.shape))
        return stamps, np.array(cols)

    def solve(dx, dy):
        # non-negative fluxes; the background may have either sign (the cutout's sky was subtracted
        # with a percentile, which in a crowded field also removes part of the faint stars' light)
        stamps, A = design(dx, dy)
        M = np.hstack([A[:, mask].T, -np.ones((mask.sum(), 1))])
        coef, rnorm = nnls(M, oot[mask])
        coef = np.append(coef[:-2], coef[-2] - coef[-1])
        return stamps, A, coef, float(rnorm ** 2)

    best = minimize(lambda p: solve(*p)[3], [0.0, 0.0], method="Nelder-Mead",
                    options={"xatol": 1e-3, "fatol": 1e-6 * np.sum(oot[mask] ** 2)})
    dx, dy = best.x
    stamps, A, coef, sse = solve(dx, dy)
    flux = np.zeros(len(xy))
    flux[free] = coef[: len(free)]
    flux[rest] = coef[len(free)] * w0[rest]
    model = np.tensordot(coef, A, axes=1)
    return {"dx": float(dx), "dy": float(dy), "flux": flux, "bkg": float(coef[-1]), "stamps": stamps,
            "model": model, "resid_frac": float(np.sqrt(sse / np.sum(oot[mask] ** 2)))}


# ------------------------------------------------------------------------------------------------
# which star is the source of the difference image?
# ------------------------------------------------------------------------------------------------


def transit_differences(t: np.ndarray, cube: np.ndarray, P: float, T0: float, dur_d: float,
                        phase_shift: float = 0.0, n_in_min: int = 3) -> np.ndarray:
    """One out-of-transit minus in-transit image per transit (local baseline on each side).

    With `phase_shift` != 0 the same is done at a fake ephemeris T0 + phase_shift * P: the
    resulting "null" difference images contain no transit and measure the noise.
    """
    T0 = T0 + phase_shift * P
    epoch = np.round((t - T0) / P).astype(int)
    dt = t - (T0 + epoch * P)
    out = []
    for k in np.unique(epoch):
        i_in = (epoch == k) & (np.abs(dt) < 0.5 * dur_d)
        i_oot = (epoch == k) & (np.abs(dt) > 0.75 * dur_d) & (np.abs(dt) < 2.5 * dur_d)
        if i_in.sum() >= n_in_min and i_oot.sum() >= 2 * n_in_min:
            out.append(np.median(cube[i_oot], axis=0) - np.mean(cube[i_in], axis=0))
    return np.array(out)


def null_differences(t, cube, P, T0, dur_d, n_transits: int, n_shift: int = 40, avoid: float = 0.06,
                     seed: int = 0, max_tries: int = 2000, n_in_min: int = 3) -> np.ndarray:
    """Mean difference images at `n_shift` fake ephemerides away from the transit (and from phase 0.5).

    Only fake ephemerides with the same number of usable transits as the real one are kept, so the
    scatter of these images is the noise of the real mean difference image.
    """
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(max_tries):
        if len(out) >= n_shift:
            break
        s = rng.uniform(avoid, 1 - avoid)
        if abs(s - 0.5) < avoid:
            continue
        d = transit_differences(t, cube, P, T0, dur_d, phase_shift=s, n_in_min=n_in_min)
        if len(d) == n_transits:
            out.append(d.mean(axis=0))
    return np.array(out)


def fit_sources(diff: np.ndarray, sigma: np.ndarray, stamps: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Weighted fit diff = a * stamp_j + b for every star j. Returns rows (a, sigma_a, chi2)."""
    w = 1.0 / sigma[mask] ** 2
    y = diff[mask]
    out = []
    for s in stamps:
        A = np.vstack([s[mask], np.ones(mask.sum())]).T
        ATA = A.T @ (A * w[:, None])
        coef = np.linalg.solve(ATA, A.T @ (w * y))
        cov = np.linalg.inv(ATA)
        chi2 = float(np.sum(w * (y - A @ coef) ** 2))
        out.append((coef[0], np.sqrt(cov[0, 0]), chi2))
    return np.array(out)


def fit_free_position(diff: np.ndarray, sigma: np.ndarray, prf: np.ndarray, x0: float, y0: float,
                      mask: np.ndarray) -> tuple[float, float, float, float]:
    """Best-fitting position of a single PRF for the difference image. Returns (x, y, a, chi2)."""
    w = 1.0 / sigma[mask] ** 2
    y = diff[mask]

    def chi2(p):
        s = prf_stamp(prf, p[0], p[1], diff.shape)[mask]
        A = np.vstack([s, np.ones_like(s)]).T
        coef = np.linalg.lstsq(A * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
        return float(np.sum(w * (y - A @ coef) ** 2)), coef[0]

    best = minimize(lambda p: chi2(p)[0], [x0, y0], method="Nelder-Mead", options={"xatol": 1e-3})
    c, a = chi2(best.x)
    return float(best.x[0]), float(best.x[1]), float(a), c


def depth_in_aperture(diff: np.ndarray, oot_model: np.ndarray, bkg: float, ap: np.ndarray) -> float:
    """Transit depth measured in aperture `ap`: lost flux over the stars' flux in the aperture."""
    return float(diff[ap].sum() / (oot_model[ap] - bkg).sum())


def predicted_depth(a: float, stamp: np.ndarray, oot_model: np.ndarray, bkg: float, ap: np.ndarray) -> float:
    """Depth an aperture would show if one star (PRF `stamp`) lost flux `a`."""
    return float(a * stamp[ap].sum() / (oot_model[ap] - bkg).sum())


def gaia_tmag(g: np.ndarray, bp_rp: np.ndarray) -> np.ndarray:
    """TESS magnitude from Gaia G and BP-RP (Stassun et al. 2019, eq. 1; T = G - 0.43 without a colour)."""
    c = np.clip(np.asarray(bp_rp, float), -0.2, 3.0)
    t = g - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473
    return np.where(np.isfinite(bp_rp), t, np.asarray(g, float) - 0.430)


def cache_dir() -> Path:
    return Path(os.environ.get("LIGHTKURVE_CACHE_DIR", "~/.lightkurve/cache")).expanduser() / "tess_prf"


# ------------------------------------------------------------------------------------------------
# the whole test for one sector
# ------------------------------------------------------------------------------------------------


def source_test(tpf, sources, P: float, T0: float, dur_d: float, prf: np.ndarray, target_id: int,
                fit_radius: float = 4.0, cand_radius: float = 3.5, n_free: int = 10, n_null: int = 60,
                seed: int = 0, epoch: float = 2000.0, n_in_min: int = 3, tmag_max_cand: float = 99.0) -> dict:
    """PRF source test on one TESScut cutout.

    `sources` is a catalogue table (columns ID, ra, dec, pmRA, pmDEC, Tmag) covering the cutout and a
    margin, with positions at Julian year `epoch` (2000 for the TIC, 2016 for Gaia DR3). Only stars
    brighter than `tmag_max_cand` are tried as the source (fainter ones stay in the scene model).
    Returns the scene fit, the mean difference image and its per-pixel noise, a table with one row per star within `cand_radius` pixels of the target (fitted lost flux, the eclipse depth that
    star would need, chi2 and delta-chi2 to the best star), and a free-position PRF fit with its
    scatter over noise realisations.
    """
    import pandas as pd
    from astropy.coordinates import SkyCoord
    from astropy.time import Time

    good = (tpf.quality == 0) & np.all(np.isfinite(tpf.flux.value), axis=(1, 2))
    t, cube = tpf.time.value[good], tpf.flux.value[good]
    cube = cube - np.percentile(cube.reshape(len(cube), -1), 20, axis=1)[:, None, None]

    src = sources.copy()
    yr = Time(np.median(t) + 2457000, format="jd").jyear - epoch
    pmra, pmde = src["pmRA"].fillna(0.0), src["pmDEC"].fillna(0.0)
    dec = src["dec"] + pmde / 3.6e6 * yr
    ra = src["ra"] + pmra / 3.6e6 * yr / np.cos(np.radians(src["dec"]))
    x, y = tpf.wcs.world_to_pixel(SkyCoord(ra.values, dec.values, unit="deg"))
    src["x"], src["y"] = x, y
    ny, nx = cube.shape[1:]
    src = src[(x > -4) & (x < nx + 3) & (y > -4) & (y < ny + 3)].reset_index(drop=True)

    dt = (t - T0 + P / 2) % P - P / 2
    oot = np.median(cube[np.abs(dt) > 1.5 * dur_d], axis=0)
    sc = scene_fit(oot, src[["x", "y"]].to_numpy(), src["Tmag"].to_numpy(), prf, n_free=n_free)
    src["flux"] = sc["flux"]
    src["x"] += sc["dx"]
    src["y"] += sc["dy"]
    it = int(src.index[src["ID"] == target_id][0])
    tx, ty = src.loc[it, "x"], src.loc[it, "y"]

    d = transit_differences(t, cube, P, T0, dur_d, n_in_min=n_in_min)
    diff = d.mean(axis=0)
    nulls = null_differences(t, cube, P, T0, dur_d, len(d), n_shift=n_null, seed=seed, n_in_min=n_in_min)
    sigma = nulls.std(axis=0)
    yy, xx = np.indices(diff.shape)
    mask = np.hypot(xx - tx, yy - ty) <= fit_radius

    near = src.index[(np.hypot(src["x"] - tx, src["y"] - ty) <= cand_radius)
                     & ((src["Tmag"] <= tmag_max_cand) | (src.index == it))]
    fits_ = fit_sources(diff, sigma, sc["stamps"][near], mask)
    tab = src.loc[near, ["ID", "Tmag", "x", "y", "flux"]].copy()
    tab["sep_px"] = np.hypot(tab["x"] - tx, tab["y"] - ty)
    tab["lost_flux"], tab["lost_flux_err"], tab["chi2"] = fits_.T
    tab["depth_needed"] = tab["lost_flux"] / tab["flux"]
    tab["dchi2"] = tab["chi2"] - tab["chi2"].min()
    tab["dchi2_vs_target"] = tab["chi2"] - tab.loc[it, "chi2"]

    # free position, and its scatter: refit (best target model + a null image) many times
    fx, fy, fa, fchi2 = fit_free_position(diff, sigma, prf, tx, ty, mask)
    tgt_model = tab.loc[it, "lost_flux"] * sc["stamps"][it]
    sim = np.array([fit_free_position(tgt_model + (n - nulls.mean(axis=0)), sigma, prf, tx, ty, mask)[:2]
                    for n in nulls])
    pos_err = sim.std(axis=0)
    return {"t": t, "cube": cube, "sources": src, "scene": sc, "diff": diff, "sigma": sigma, "nulls": nulls,
            "n_transits": len(d), "mask": mask, "target": (tx, ty), "table": tab.sort_values("chi2"),
            "free": {"x": fx, "y": fy, "lost_flux": fa, "chi2": fchi2, "x_err": float(pos_err[0]),
                     "y_err": float(pos_err[1]),
                     "offset_px": float(np.hypot(fx - tx, fy - ty)),
                     "offset_sigma": float(np.hypot((fx - tx) / pos_err[0], (fy - ty) / pos_err[1]))},
            "n_pix": int(mask.sum())}
