"""End-to-end on LHS 3844 (TIC 410153553) with the real sector-104 data: search must find
LHS 3844 b, vetting must pass it, and the known-object step must label it a known planet.

Needs data/ from a full run (skipped otherwise) and runs TLS once (~1 min): `pytest -m slow`.
"""

from __future__ import annotations

import os
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

TIC = 410153553
P_TRUE = 0.462930  # Vanderspek et al. 2019
DATA = Path(os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data"))
SECTOR = int(os.environ.get("SECTOR", 104))

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not (DATA / "lc_flat" / f"{TIC}.parquet").exists(), reason="sector data not on disk"),
]


@pytest.fixture(scope="module")
def search_row():
    from src.search import search_one

    targets = pd.read_parquet(DATA / "targets" / f"sector_{SECTOR}.parquet").set_index("tic_id")
    row = search_one(TIC, str(DATA / "lc_flat" / f"{TIC}.parquet"), targets.loc[TIC, "rad"], targets.loc[TIC, "mass"])
    assert row.status == "ok", row.error
    return row


def test_search_finds_lhs3844b(search_row):
    assert abs(search_row.period - P_TRUE) / P_TRUE < 1e-3
    assert search_row.SDE > 20


def test_vetting_passes_lhs3844b(search_row):
    from src.vet import apply_thresholds, momentum_dumps, vet_one

    targets = pd.read_parquet(DATA / "targets" / f"sector_{SECTOR}.parquet").set_index("tic_id")
    manifest = pd.read_parquet(DATA / "targets" / f"sector_{SECTOR}_manifest.parquet")
    s = asdict(search_row)
    v = vet_one(TIC, s, targets.loc[TIC].to_dict(), str(DATA / "lc_flat" / f"{TIC}.parquet"),
                manifest.set_index("tic_id").loc[TIC, "file_path"], momentum_dumps(manifest))
    row = {**asdict(v), "SDE": s["SDE"], "snr_tls": s["snr"], "distinct_transit_count": s["distinct_transit_count"],
           "period": s["period"], "period_max": s["period_max"], "period_unc": s["period_unc"]}
    out = apply_thresholds(pd.DataFrame([row])).iloc[0]
    assert out["passes_all"], out["failed"]
    assert abs(out["rp_earth"] - 1.3) < 0.3  # published 1.30 R_earth


def test_known_step_labels_it_known_planet():
    path = DATA / "search" / f"sector_{SECTOR}_known.parquet"
    if not path.exists():
        pytest.skip("run src.known first")
    k = pd.read_parquet(path).set_index("tic_id").loc[TIC]
    assert k["known_class"] == "known_planet"
    assert "LHS 3844 b" in k["known_ids"]
    assert np.isfinite(k["snr"])
