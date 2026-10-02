"""Matching math of src.known, including regression tests for bugs found during development."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd

from src import known
from src.known import depth_needed_ppm, ephemeris_check, likely_source, period_ratio


def cat_row(**kw):
    base = dict(catalog="planet", period=2.0, period_err=1e-6, epoch=1500.0, epoch_err=1e-4)
    return SimpleNamespace(**{**base, **kw})


def test_period_ratio_aliases():
    assert period_ratio(2.0, 2.001)[0] == 1
    assert period_ratio(1.0, 2.0)[0] == 0.5
    assert period_ratio(6.0, 2.0)[0] == 3
    assert period_ratio(2.0, 2.1) is None  # 5 % off: no match
    assert period_ratio(2.0, np.nan) is None


def test_ephemeris_uses_catalog_period():
    """Regression (2026-10-01): propagating with OUR period over ~1,300 orbits turned a 0.002-d
    period difference into hours, and three confirmed planets were called new."""
    k = cat_row(period=2.0, epoch=1500.0)
    n = 1340
    true_T0 = 1500.0 + n * 2.0
    dt, sigma_t, status = ephemeris_check(true_T0 + 0.003, 1.0, k, 1)  # our epoch, 4 minutes late
    assert status == "ok"
    assert abs(dt - 0.003) < 1e-9
    assert abs(sigma_t - np.hypot(1e-4, n * 1e-6)) < 1e-9


def test_ephemeris_rejects_wrong_phase_and_handles_aliases():
    k = cat_row(period=2.0, epoch=1500.0)
    assert ephemeris_check(1500.0 + 100 * 2.0 + 0.7, 1.0, k, 1)[2] == "bad"
    # our period = half the catalog's (EB folded at P/2): dips also at the catalog's phase 0.5
    assert ephemeris_check(1500.0 + 100 * 2.0 + 1.0, 1.0, k, 0.5)[2] == "ok"
    # no epoch error -> unknown, never "bad"
    assert ephemeris_check(1600.7, 1.0, cat_row(epoch_err=np.nan), 1)[2] == "unknown"


def test_eb_secondary_counts_as_match():
    k = cat_row(catalog="tess_eb", period=2.0, epoch=1500.0)
    assert ephemeris_check(1500.0 + 50 * 2.0 + 1.0, 1.0, k, 1)[2] == "ok"


def test_dilution_deblended_formula():
    """A neighbour 2.5 mag brighter (10x the flux) needs only a tenth of our depth."""
    assert abs(depth_needed_ppm(5000, 12.5, 10.0) - 500) < 1e-6
    assert likely_source(5000, 12.5, 10.0, 600)[1] == "neighbour"
    assert likely_source(5000, 12.5, 10.0, 100)[1] == "target"
    assert likely_source(5000, 10.0, 17.5, 1e5)[1] == "target"  # neighbour far too faint
    assert likely_source(5000, 12.5, 10.0, np.nan)[1] == "unclear"


def test_gaia_offline_merges_caches_of_all_days(tmp_path, monkeypatch):
    """Regression (2026-10-02): offline mode read only the newest dated folder, found no Gaia
    cache there and silently skipped the Gaia check (8 blends became 'new')."""
    monkeypatch.setattr(known, "KNOWN_ROOT", tmp_path)
    old, new = tmp_path / "2026-10-01", tmp_path / "2026-10-02"
    old.mkdir()
    new.mkdir()
    pd.DataFrame({"tic_query": [11], "Source": [123], "RA_ICRS": [10.0], "DE_ICRS": [-40.0], "Freq": [0.5]}).to_parquet(
        old / "gaia_veb_sector_1.parquet", index=False)
    (old / "gaia_veb_sector_1_queried.json").write_text(json.dumps([11, 12]))
    sig = pd.DataFrame({"ra": [10.01, 50.0], "dec": [-40.0, 10.0], "gaia_dr2_id": [999, 998]}, index=[11, 12])
    g = known.gaia_ebs(sig, new, 1, offline=True)
    assert len(g) == 1 and g["tic_query"].iloc[0] == 11
    assert abs(g["period"].iloc[0] - 2.0) < 1e-12
    assert not g["self"].iloc[0]  # 36" away and a different Gaia ID
