"""Follow-up repetition check: the narrow-window search must stay narrow; weak sectors cannot refute."""

from __future__ import annotations

import numpy as np
import pandas as pd

from conftest import CADENCE, MASS, NOISE, RAD, time_grid
from src.followup import SNR_INFORMATIVE, _tls_window, expected_snr, repetition
from src.inject import transit_model


def test_tls_window_stays_in_window(rng):
    """One sector of a 4.5-d planet: the +-2 % window holds < 100 default grid periods, which made TLS fall
    back to a full search. The window must be respected and the planet recovered."""
    t = time_grid()
    P, T0 = 4.5, 4181.0
    model, depth, dur_h = transit_model(t, P, T0, 3.0, 0.2, RAD, MASS, CADENCE)
    lc = pd.DataFrame({"time": t, "flux": model + rng.normal(0, NOISE, len(t)), "flux_err": NOISE})
    r, ovs = _tls_window(lc, P, RAD, MASS, 0.02)
    assert ovs > 3
    assert abs(r.period - P) / P < 0.005
    assert expected_snr(lc, P, dur_h, depth * 1e6) > SNR_INFORMATIVE


def _sec_table(rows):
    return pd.DataFrame(rows, columns=["data", "period", "snr", "depth_ppm", "snr_expected"])


def test_weak_sectors_cannot_refute():
    P = 2.0
    weak = _sec_table([("s104", P, 8.0, 3000, 8.0), ("s001", 1.3, 3.0, 2000, 3.0), ("s028", 0.7, 3.5, 2500, 4.0)])
    r = repetition(weak, P, 3000, 104)
    assert r["n_informative"] == 0 and not r["testable"] and not r["repeats"]
    # ... but the weak sectors together are informative: the combined search decides
    absent = pd.concat([weak, _sec_table([("all but discovery", 2.07, 0.1, 900, 12.0)])], ignore_index=True)
    r = repetition(absent, P, 3000, 104)
    assert r["basis"] == "combined" and r["testable"] and not r["repeats"]
    present = pd.concat([weak, _sec_table([("all but discovery", P * 1.0003, 0.1, 2700, 12.0)])], ignore_index=True)
    assert repetition(present, P, 3000, 104)["repeats"]
    strong = _sec_table([("s104", P, 8.0, 3000, 8.0), ("s101", 1.3, 3.0, 2000, 9.0), ("s102", 0.7, 3.5, 2500, 9.0)])
    r = repetition(strong, P, 3000, 104)
    assert r["n_informative"] == 2 and r["testable"] and not r["repeats"]
    real = _sec_table([("s104", P, 8.0, 3000, 8.0), ("s101", P * 1.001, 7.0, 2800, 9.0),
                       ("s102", P * 0.999, 6.0, 3100, 9.0), ("s001", 1.1, 3.0, 2000, 3.0)])
    assert repetition(real, P, 3000, 104)["repeats"]


def test_eb_at_half_period_is_deeper_at_2P():
    P = 0.6
    rows = [("s104", P, 9.0, 4300, 9.0), ("s066", P, 16.0, 11000, 9.0), ("s093", P, 12.0, 7200, 9.0),
            ("s066@2P", 2 * P, 14.0, 14000, np.nan), ("s093@2P", 2 * P, 10.0, 9100, np.nan)]
    assert repetition(_sec_table(rows), P, 4300, 104)["repeats_at_alias"]
    planet = [("s104", P, 9.0, 4300, 9.0), ("s066", P, 10.0, 4200, 9.0), ("s093", P, 9.0, 4400, 9.0),
              ("s066@2P", 2 * P, 7.0, 4300, np.nan), ("s093@2P", 2 * P, 6.5, 4250, np.nan)]
    r = repetition(_sec_table(planet), P, 4300, 104)
    assert r["repeats"] and not r["repeats_at_alias"]
