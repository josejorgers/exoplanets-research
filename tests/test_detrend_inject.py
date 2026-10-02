"""Detrending must preserve transits; the injection model must have the right depth and duration."""

from __future__ import annotations

import numpy as np

from conftest import CADENCE, MASS, RAD, time_grid
from src.detrend import edge_mask, flatten_lc
from src.inject import transit_model
from src.vet import _expected_duration_h


def test_edge_mask_trims_each_orbit():
    t = time_grid()
    keep = edge_mask(t, 0.25)
    for lo, hi in [(t.min(), 4192.0), (4193.0, t.max())]:
        seg = (t >= lo) & (t <= hi)
        kept = t[seg & keep]
        assert kept.min() - t[seg].min() >= 0.25 - 1e-9
        assert t[seg].max() - kept.max() >= 0.25 - 1e-9


def test_flatten_preserves_a_2mmag_box_on_a_trend():
    """Plan §5 sanity check: a 2-mmag, 1-h box on a 1 %, 3-d sine trend with 1,000 ppm noise keeps
    >= 80 % of its depth through the 0.5-d biweight (measured 0.90; 2.5-h boxes keep only ~0.75,
    see docs/NOTES_decisions.md). Noiseless data would be a misleading test: the biweight's MAD
    collapses and the trend's curvature biases it."""
    t = time_grid()
    rng = np.random.default_rng(0)
    box = np.where(np.abs((t - 4181.3 + 1.0) % 2.0 - 1.0) < 0.5 / 24, 1 - 0.002, 1.0)
    flux = (1 + 0.01 * np.sin(2 * np.pi * t / 3.0)) * box + rng.normal(0, 1e-3, len(t))
    out, _ = flatten_lc(t, flux, np.full_like(t, 1e-3))
    tt, y = out["time"].to_numpy(), out["flux"].to_numpy()
    in_box = np.abs((tt - 4181.3 + 1.0) % 2.0 - 1.0) < 0.33 / 24
    retained = (1 - np.median(y[in_box]) / np.median(y[~in_box])) / 0.002
    assert 0.8 < retained < 1.1
    assert 1.4826 * np.median(np.abs(y[~in_box] - 1)) < 1.2e-3  # the 1 % trend is gone: only the noise is left


def test_transit_model_depth_and_duration():
    t = np.arange(4180.0, 4182.0, CADENCE / 10)
    k = 3.0 / (RAD * 109.08)
    flux, depth, dur_h = transit_model(t, 1.5, 4181.0, 3.0, 0.0, RAD, MASS, CADENCE / 10)
    assert 1.0 * k**2 < depth < 1.4 * k**2  # limb darkening deepens a central transit
    assert abs(dur_h - _expected_duration_h(1.5, RAD, MASS, k)) / dur_h < 0.02
    assert flux.min() > 1 - depth - 1e-4
