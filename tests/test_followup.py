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


def test_trilegal_roundtrip_for_triceratops(tmp_path):
    """Our TRILEGAL parser + writer must give a file triceratops' own reader accepts, with no star lost
    (it drops the last two rows, which its own writer fills with the termination comment)."""
    from triceratops.funcs import trilegal_results

    from src.background import parse_trilegal, write_triceratops_csv

    cols = "#Gc logAge [M/H] m_ini logL logTe logg m-M0 Av m2/m1 mbol TESS J H Ks Mact"
    rows = ["1 9.1 0.0 0.40 -1.5 3.55 4.8 10.2 0.1 0.0 14.6 15.5 14.0 13.4 13.2 0.40",
            "1 9.5 -0.2 0.80 -0.3 3.70 4.6 12.0 0.2 0.0 13.5 16.8 15.6 15.2 15.1 0.80",
            "2 10.0 -0.5 0.90 0.1 3.72 4.4 14.0 0.3 0.0 14.0 18.9 17.8 17.4 17.3 0.89"]
    text = "\n".join([cols, *rows, "#TRILEGAL normally terminated"]) + "\n"
    df = parse_trilegal(text)
    assert len(df) == 3 and df["TESS"].tolist() == [15.5, 16.8, 18.9]
    f = tmp_path / "trilegal.csv"
    write_triceratops_csv(df, f)
    tmags = trilegal_results(str(f), 16.0)[0]
    assert np.allclose(tmags, [16.8, 18.9])          # only stars fainter than the target, none dropped


def test_trilegal_incomplete_output_rejected():
    import pytest

    from src.background import parse_trilegal

    with pytest.raises(RuntimeError, match="incomplete"):
        parse_trilegal("#Gc TESS J H Ks Mact logg logTe [M/H]\n1 15 14 13 13 0.4 4.8 3.55 0.0\n")
    with pytest.raises(RuntimeError, match="lacks columns"):
        parse_trilegal("#Gc J H Ks\n1 14 13 13\n#TRILEGAL normally terminated\n")


def test_scenario_groups_partition_triceratops_scenarios():
    """Every TRICERATOPS scenario belongs to exactly one origin group, so the groups sum to 1."""
    from src.followup import SCENARIO_GROUPS, scenario_groups

    names = ["TP", "EB", "EBx2P", "PTP", "PEB", "PEBx2P", "STP", "SEB", "SEBx2P", "DTP", "DEB", "DEBx2P",
             "BTP", "BEB", "BEBx2P", "NTP", "NEB", "NEBx2P"]
    flat = [s for sc in SCENARIO_GROUPS.values() for s in sc]
    assert sorted(flat) == sorted(set(names))
    probs = pd.DataFrame({"scenario": names, "ID": "1", "prob": np.full(len(names), 1 / len(names))})
    g = scenario_groups(probs)
    assert abs(sum(g.values()) - 1) < 1e-12
    assert abs(g["background"] - 3 / len(names)) < 1e-12
