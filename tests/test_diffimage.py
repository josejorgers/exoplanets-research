"""Difference-image decision on synthetic cutouts: on target, off target, and noise with few transits."""

import numpy as np

from src.followup import difference_image_stats

SIZE, TX, TY = 15, 7.0, 7.0
P, T0, DUR = 3.0, 1.0, 2.0 / 24
CADENCE = 200 / 86400


def _cube(n_days, depth=0.0, at=(TX, TY), noise=30.0, seed=1, flux=2.0e4, hot=None):
    """A Gaussian star at the target (and a neighbour 4 px away) plus white noise; `depth` dims the star at `at`."""
    rng = np.random.default_rng(seed)
    t = np.arange(0, n_days, CADENCE)
    yy, xx = np.indices((SIZE, SIZE))

    def psf(x0, y0):
        g = np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * 0.9 ** 2))
        return g / g.sum()

    target, neighbour = flux * psf(TX, TY), 0.5 * flux * psf(TX + 4, TY)
    in_tr = np.abs((t - T0 + 0.5 * P) % P - 0.5 * P) < 0.5 * DUR
    cube = np.broadcast_to(target + neighbour, (len(t), SIZE, SIZE)).copy()
    dimmed = psf(*at) * (flux if at == (TX, TY) else 0.5 * flux)
    cube[in_tr] -= depth * dimmed
    cube += rng.normal(0, noise, cube.shape)
    if hot is not None:  # one pixel with extra noise, the kind that wins a per-pixel S/N map
        cube[:, hot[1], hot[0]] += rng.normal(0, 4 * noise, len(t))
    return t, cube


def test_transit_on_target():
    t, cube = _cube(27, depth=0.01)
    r = difference_image_stats(t, cube, P, T0, DUR, TX, TY)
    assert r["status"] == "on_target" and not r["off_target"]
    assert r["noise_source"] == "pseudo-transits"
    assert r["offset_from_catalog_px"] < 0.3


def test_transit_on_neighbour():
    t, cube = _cube(27, depth=0.02, at=(TX + 4, TY))
    r = difference_image_stats(t, cube, P, T0, DUR, TX, TY)
    assert r["status"] == "off_target" and r["off_target"]
    assert abs(r["strongest_loss_dist_px"] - 4) < 1


def test_noise_with_few_transits_is_not_off_target():
    # 4 transits, no signal, one noisy pixel far from the target: nothing should be called off target
    for seed in range(5):
        t, cube = _cube(12, depth=0.0, seed=seed, hot=(13, 2))
        r = difference_image_stats(t, cube, P, T0, DUR, TX, TY)
        assert r["n_transits"] == 4
        assert r["status"] == "inconclusive" and not r["off_target"], (seed, r["strongest_loss_snr"], r["snr_threshold"])


def test_weak_transit_few_events_still_on_target():
    t, cube = _cube(12, depth=0.004)
    r = difference_image_stats(t, cube, P, T0, DUR, TX, TY)
    assert r["status"] == "on_target"


def test_too_few_transits():
    t, cube = _cube(4.5, depth=0.01)
    r = difference_image_stats(t, cube, 5.0, 4.0, DUR, TX, TY)
    assert "error" in r
