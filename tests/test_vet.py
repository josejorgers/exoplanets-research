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


def test_two_transit_planet_passes_single_event(tmp_path, star, rng):
    """Regression (2026-10-02, found by injection-recovery): with 2 transits, dropping the stronger
    one left fewer than 2 and the test returned 0, rejecting every planet with P > ~8 d."""
    t = time_grid()
    P, T0 = 11.0, 4181.5
    model, depth, dur_h = transit_model(t, P, T0, 3.0, 0.2, RAD, MASS, CADENCE)
    flux = model + rng.normal(0, NOISE, len(t))
    sig = make_sig(P, T0, dur_h, depth, 3.0 / (RAD * 109.08), n_transits=2)
    out = vet_table(sig, vet_one(1, sig, star, write_lc(tmp_path, t, flux), None, NO_DUMPS)).iloc[0]
    assert out["n_transits_used"] == 2
    assert out["single_event_survival"] > 0.5
    assert out["pass_single_event"], out["failed"]


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


def test_transits_at_data_edges_fail(tmp_path, star, rng):
    """'Transits' that sit right at the start of data segments (ramps after gaps) are not trusted."""
    t = time_grid()
    t = t[((t - 4180.0) % 3.0) > 0.4]  # 0.4-d gaps every 3 days
    flux = 1 + rng.normal(0, NOISE, len(t))
    starts = t[np.r_[0, np.where(np.diff(t) > 0.1)[0] + 1]]
    for s0 in starts:  # a dip in the first hour of every segment
        flux[(t >= s0) & (t < s0 + 1 / 24)] -= 0.004
    sig = make_sig(3.0, starts[1] + 0.5 / 24, 1.0, 0.004, 0.06, n_transits=8)
    out = vet_table(sig, vet_one(1, sig, star, write_lc(tmp_path, t, flux), None, NO_DUMPS)).iloc[0]
    assert out["n_edge"] >= 6
    assert not out["pass_data_edges"]


def test_noise_bursts_inside_transit_windows_fail(tmp_path, star, rng):
    """Bursts of 8x scatter confined to 3 of the predicted transit windows (like TIC 394134586):
    their mean is low, so TLS sees dips, but the in-transit scatter gives them away."""
    t = time_grid()
    flux = 1 + rng.normal(0, NOISE, len(t))
    P, T0, dur_h = 2.2, 4180.9, 1.8
    for n in (8, 9, 10):
        m = np.abs(t - (T0 + n * P)) < 0.6 / 24
        flux[m] = 1 - 0.006 + rng.normal(0, 8 * NOISE, m.sum())
    sig = make_sig(P, T0, dur_h, 0.0015, 0.04, n_transits=10)
    out = vet_table(sig, vet_one(1, sig, star, write_lc(tmp_path, t, flux), None, NO_DUMPS)).iloc[0]
    assert out["n_noisy"] >= 3
    assert not out["pass_transit_noise"]


def test_local_noise_downweights_transit_in_noisy_stretch(planet, star, tmp_path):
    """The same planet, with one day of 10x noise around one transit: that transit's error grows."""
    sig, path, (P, T0, dur_h, depth) = planet
    lc = pd.read_parquet(path)
    noisy = np.abs(lc["time"] - (T0 + 3 * P)) < 0.5
    rng = np.random.default_rng(1)
    lc.loc[noisy, "flux"] += rng.normal(0, 10 * NOISE, noisy.sum())
    lc.to_parquet(tmp_path / "noisy.parquet", index=False)
    from src.vet import _per_transit, _phase_time

    t, y = lc["time"].to_numpy(), lc["flux"].to_numpy()
    dt = _phase_time(t, P, T0)
    dur = dur_h / 24
    pt = _per_transit(t, y, np.abs(dt) < dur / 2, np.abs(dt) > dur, np.round((t - T0) / P).astype(int),
                      T0, P, dur, 1.0, NOISE)
    err = pt.set_index("epoch")["err"]
    assert err.loc[3] > 5 * err.drop(3).median()


def test_variability_matches_harmonics_in_frequency(tmp_path, star, rng):
    """A 2.10-d spot modulation found by TLS at 2.09 d (0.5 % off) is still recognised: the
    periodogram only locates the peak to ~1/baseline in frequency (~8 % in period here)."""
    t = time_grid()
    flux = 1 + 0.003 * np.sin(2 * np.pi * t / 2.10) + rng.normal(0, NOISE, len(t))
    sig = make_sig(2.09, 4180.0 + 0.75 * 2.10, 2.4, 0.004, 0.05)
    out = vet_table(sig, vet_one(1, sig, star, write_lc(tmp_path, t, flux), None, NO_DUMPS)).iloc[0]
    assert out["var_fap"] < 1e-3
    assert out["var_harmonic_offset"] < 1
    assert not out["pass_variability"]


def test_planet_on_quiet_star_is_not_called_variable(planet, star):
    sig, path, _ = planet
    out = vet_table(sig, vet_one(1, sig, star, path, None, NO_DUMPS)).iloc[0]
    assert out["var_fap"] > 1e-3  # nothing but white noise once the transits are masked
    assert out["pass_variability"] and not out["flag_variable_star"]


def test_nan_metric_passes_by_convention():
    """A test that could not be evaluated must not reject the signal."""
    row = {"status": "ok", "SDE": 20.0, "snr_tls": 20.0, "distinct_transit_count": 5, "period": 1.0,
           "period_max": 12.0, "period_unc": 0.001, "snr": 20.0}
    metrics = ["single_event_survival", "dump_survival", "edge_survival", "noisy_survival", "odd_even_sigma", "secondary_snr",
               "duration_ratio", "rp_earth", "var_fap", "var_harmonic_offset", "var_amp_frac", "trap_flat_frac",
               "centroid_sigma", "contratio"]
    out = apply_thresholds(pd.DataFrame([{**row, **{m: np.nan for m in metrics}}])).iloc[0]
    assert out["passes_all"]
    assert out["n_flags"] == 0


def test_expected_duration_earth_sun():
    """An Earth-Sun central transit lasts ~13 hours."""
    assert abs(_expected_duration_h(365.25, 1.0, 1.0) - 13.0) < 0.3
