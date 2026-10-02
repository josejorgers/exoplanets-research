"""Each vetting test must pass a clean synthetic planet and reject the impostor it is designed for."""

from __future__ import annotations

import numpy as np
import pandas as pd

from conftest import NOISE, RAD, MASS, CADENCE, make_sig, time_grid, write_lc
from src.inject import transit_model
from src.vet import _expected_duration_h, apply_thresholds, vet_one

NO_DUMPS = np.array([])


def vet_table(sig: dict, row) -> pd.DataFrame:
    """One vetted row in the format apply_thresholds expects (search columns + metrics)."""
    from dataclasses import asdict

    d = {**asdict(row), "SDE": sig["SDE"], "snr_tls": 30.0, "distinct_transit_count": sig["distinct_transit_count"],
         "period": sig["period"], "period_max": 12.0, "period_unc": 0.001}
    return apply_thresholds(pd.DataFrame([d]))


def test_planet_passes_everything(planet, star):
    sig, path, _ = planet
    row = vet_one(1, sig, star, path, None, NO_DUMPS)
    assert row.status == "ok", row.error
    out = vet_table(sig, row).iloc[0]
    assert out["passes_all"], out["failed"]
    assert out["odd_even_sigma"] < 3
    assert out["secondary_snr"] < 5
    assert 0.5 < out["duration_ratio"] < 1.5
    assert 1.5 < out["rp_earth"] < 3.5


def test_eb_at_half_period_fails_odd_even(tmp_path, star, rng):
    """Two eclipses per orbit of different depth, folded at half the true period."""
    t = time_grid()
    P_true, T0 = 3.0, 4180.5
    prim, d1, dur = transit_model(t, P_true, T0, 4.0, 0.2, RAD, MASS, CADENCE)
    sec, _, _ = transit_model(t, P_true, T0 + P_true / 2, 2.5, 0.2, RAD, MASS, CADENCE)
    flux = prim * sec + rng.normal(0, NOISE, len(t))
    sig = make_sig(P_true / 2, T0, dur, d1 * 0.7, 0.1)
    row = vet_one(1, sig, star, write_lc(tmp_path, t, flux), None, NO_DUMPS)
    out = vet_table(sig, row).iloc[0]
    assert out["odd_even_sigma"] > 3
    assert not out["pass_odd_even"]


def test_secondary_eclipse_fails(tmp_path, star, rng):
    t = time_grid()
    P, T0 = 2.0, 4180.4
    prim, depth, dur = transit_model(t, P, T0, 4.0, 0.2, RAD, MASS, CADENCE)
    dt = (t - T0 - P / 2 + P / 2) % P - P / 2
    secondary = np.where(np.abs(dt) < dur / 48, 1 - 0.3 * depth, 1.0)  # 30 % of the primary, at phase 0.5
    flux = prim * secondary + rng.normal(0, NOISE, len(t))
    sig = make_sig(P, T0, dur, depth, 0.12)
    out = vet_table(sig, vet_one(1, sig, star, write_lc(tmp_path, t, flux), None, NO_DUMPS)).iloc[0]
    assert out["secondary_snr"] > 5
    assert abs(out["secondary_phase"] - 0.5) < 0.05
    assert not out["pass_secondary"]


def test_sinusoid_fails_variability(tmp_path, star, rng):
    """Spot modulation with P_var = 0.4 d; TLS 'finds' a dip every 2 cycles."""
    t = time_grid()
    p_var, amp = 0.4, 3000e-6
    flux = 1 + amp * np.sin(2 * np.pi * t / p_var) + rng.normal(0, NOISE, len(t))
    T0 = 4180.0 + 0.75 * p_var  # a minimum of the sine
    sig = make_sig(2 * p_var, T0, 2.0, 2 * amp, 0.05)
    out = vet_table(sig, vet_one(1, sig, star, write_lc(tmp_path, t, flux), None, NO_DUMPS)).iloc[0]
    assert abs(out["var_period"] - p_var) / p_var < 0.01
    assert not out["pass_variability"]


def test_single_event_fails(tmp_path, star, rng):
    """One deep 2-hour dip; the 'period' puts the other transits on noise."""
    t = time_grid()
    flux = 1 + rng.normal(0, NOISE, len(t))
    flux[np.abs(t - 4185.3) < 1 / 24] -= 0.01
    sig = make_sig(3.0, 4185.3 - 3.0, 2.0, 0.002, 0.05, n_transits=8)
    out = vet_table(sig, vet_one(1, sig, star, write_lc(tmp_path, t, flux), None, NO_DUMPS)).iloc[0]
    assert out["max_ses_frac"] > 0.8
    assert not out["pass_single_event"]


def test_momentum_dump_dominated_signal_fails(tmp_path, star, rng):
    """'Transits' that exist only next to momentum dumps."""
    t = time_grid()
    dumps = np.array([4181.0, 4187.0, 4196.0, 4202.0])
    flux = 1 + rng.normal(0, NOISE, len(t))
    for d in dumps:
        flux[np.abs(t - d) < 1 / 24] -= 0.004
    sig = make_sig(3.0, 4181.0, 2.0, 0.004 / 3, 0.04, n_transits=8)
    out = vet_table(sig, vet_one(1, sig, star, write_lc(tmp_path, t, flux), None, dumps)).iloc[0]
    assert out["n_near_dump"] >= 2
    assert not out["pass_momentum_dump"]


def test_nan_metric_passes_by_convention():
    """A test that could not be evaluated must not reject the signal."""
    row = {"status": "ok", "SDE": 20.0, "snr_tls": 20.0, "distinct_transit_count": 5, "period": 1.0,
           "period_max": 12.0, "period_unc": 0.001, "snr": 20.0}
    metrics = ["single_event_survival", "dump_survival", "odd_even_sigma", "secondary_snr", "duration_ratio",
               "rp_earth", "var_period_ratio", "var_amp_frac", "trap_flat_frac", "centroid_sigma", "contratio"]
    out = apply_thresholds(pd.DataFrame([{**row, **{m: np.nan for m in metrics}}])).iloc[0]
    assert out["passes_all"]
    assert out["n_flags"] == 0


def test_expected_duration_earth_sun():
    """An Earth-Sun central transit lasts ~13 hours."""
    assert abs(_expected_duration_h(365.25, 1.0, 1.0) - 13.0) < 0.3
