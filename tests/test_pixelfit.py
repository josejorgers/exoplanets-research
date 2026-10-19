"""PRF pixel test (src/pixelfit.py) on synthetic images with a Gaussian stand-in for the TESS PRF."""

import numpy as np

from src import pixelfit as pf

SHAPE = (11, 11)


def gaussian_prf(sigma_px: float = 0.8) -> np.ndarray:
    n = 13 * pf.PRF_SUB
    c = (n - 1) / 2
    yy, xx = np.indices((n, n))
    return np.exp(-0.5 * (np.hypot(xx - c, yy - c) / (sigma_px * pf.PRF_SUB)) ** 2)


def test_prf_stamp_normalised_and_centred():
    prf = gaussian_prf()
    s = pf.prf_stamp(prf, 5.3, 4.6, SHAPE)
    assert abs(s.sum() - 1) < 0.02
    yy, xx = np.indices(SHAPE)
    assert abs((s * xx).sum() / s.sum() - 5.3) < 0.05
    assert abs((s * yy).sum() / s.sum() - 4.6) < 0.05
    # a star near the edge keeps only part of its light in the cutout
    assert pf.prf_stamp(prf, 0.0, 5.0, SHAPE).sum() < 0.75


def test_scene_fit_recovers_fluxes_and_offset():
    prf = gaussian_prf()
    xy = np.array([[5.0, 5.0], [6.2, 5.1], [2.0, 8.0]])
    flux = np.array([1000.0, 600.0, 300.0])
    img = sum(f * pf.prf_stamp(prf, x + 0.2, y - 0.1, SHAPE) for f, (x, y) in zip(flux, xy)) + 5.0
    sc = pf.scene_fit(img, xy, np.array([13.0, 13.5, 14.3]), prf, n_free=3)
    assert abs(sc["dx"] - 0.2) < 0.03 and abs(sc["dy"] + 0.1) < 0.03
    assert np.allclose(sc["flux"], flux, rtol=0.03)
    assert sc["resid_frac"] < 0.01


def test_fit_sources_prefers_true_source_and_excludes_neighbour():
    prf = gaussian_prf()
    rng = np.random.default_rng(3)
    xy = [(5.0, 5.0), (5.55, 5.15), (5.1, 5.0)]  # target, a neighbour 0.57 px away, one 0.1 px away
    stamps = np.array([pf.prf_stamp(prf, x, y, SHAPE) for x, y in xy])
    sigma = np.full(SHAPE, 0.5)
    diff = 15.0 * stamps[0] + rng.normal(0, 0.5, SHAPE)
    yy, xx = np.indices(SHAPE)
    mask = np.hypot(xx - 5, yy - 5) <= 4
    a, a_err, chi2 = pf.fit_sources(diff, sigma, stamps, mask).T
    assert abs(a[0] - 15) < 4 * a_err[0]
    assert chi2[1] - chi2[0] > 9  # the 0.57-px neighbour is excluded at > 3 sigma
    assert abs(chi2[2] - chi2[0]) < chi2[1] - chi2[0]  # a 0.1-px neighbour is much harder to exclude
    x, y, _, _ = pf.fit_free_position(diff, sigma, prf, 5.3, 5.3, mask)
    assert np.hypot(x - 5, y - 5) < 0.15


def test_transit_differences_and_nulls():
    prf = gaussian_prf()
    stamp = pf.prf_stamp(prf, 5.0, 5.0, SHAPE)
    t = np.arange(0, 20, 0.01)
    P, T0, dur = 4.0, 1.0, 0.1
    dt = (t - T0 + P / 2) % P - P / 2
    rng = np.random.default_rng(0)
    cube = 1000.0 * stamp[None] * (1 - 0.01 * (np.abs(dt) < dur / 2))[:, None, None]
    cube = cube + rng.normal(0, 0.05, cube.shape)
    d = pf.transit_differences(t, cube, P, T0, dur)
    assert len(d) == 5
    assert abs(d.mean(axis=0).sum() - 10.0) < 0.5  # 1 % of 1000
    nulls = pf.null_differences(t, cube, P, T0, dur, len(d), n_shift=10)
    assert len(nulls) == 10
    assert abs(nulls.mean(axis=0).sum()) < 1.0


def test_gaia_tmag():
    # Stassun et al. 2019: T is ~0.4-1.2 mag brighter than G for K and M dwarfs; G - 0.43 without a colour
    t = pf.gaia_tmag(np.array([14.28, 14.0]), np.array([2.756, np.nan]))
    assert 12.9 < t[0] < 13.3
    assert abs(t[1] - 13.57) < 1e-6
