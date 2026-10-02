"""Shared fixtures: synthetic flattened light curves with known signals.

Each fixture returns (sig, star, path): `sig` is the dict a search row would provide, `star` the
target row, `path` a parquet file in the same format as data/lc_flat/{tic}.parquet.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.inject import transit_model

RAD, MASS = 0.30, 0.30  # a typical M dwarf, solar units
NOISE = 1000e-6  # 1,000 ppm per 2-min point
CADENCE = 2 / 1440


def time_grid() -> np.ndarray:
    """25 days of 2-min cadence with a 1-day mid-sector gap, like a TESS sector."""
    t = np.arange(4180.0, 4205.0, CADENCE)
    return t[(t < 4192.0) | (t > 4193.0)]


def write_lc(tmp_path, t, flux, name="lc") -> str:
    path = tmp_path / f"{name}.parquet"
    pd.DataFrame({"time": t, "flux": flux, "flux_err": np.full_like(t, NOISE), "trend": 1.0,
                  "flux_raw": flux}).to_parquet(path, index=False)
    return str(path)


def make_sig(period, T0, dur_h, depth, rp_rs, n_transits=10, sde=20.0):
    return {"period": period, "T0": T0, "duration_h": dur_h, "depth_ppm": depth * 1e6, "rp_rs": rp_rs,
            "SDE": sde, "distinct_transit_count": n_transits}


@pytest.fixture
def star():
    return {"rad": RAD, "mass": MASS, "contratio": 0.01}


@pytest.fixture
def rng():
    return np.random.default_rng(42)


@pytest.fixture
def planet(tmp_path, rng):
    """A 2.5 R_earth planet on a 1.7-d orbit: every vetting test must pass it."""
    t = time_grid()
    P, T0, rp = 1.7, 4180.6, 2.5
    model, depth, dur_h = transit_model(t, P, T0, rp, 0.3, RAD, MASS, CADENCE)
    flux = model + rng.normal(0, NOISE, len(t))
    k = rp / (RAD * 109.08)
    return make_sig(P, T0, dur_h, depth, k), write_lc(tmp_path, t, flux, "planet"), (P, T0, dur_h, depth)
