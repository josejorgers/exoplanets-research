"""Neighbour parameters for TRICERATOPS (src/neighbours.py): no network."""

import json

import numpy as np
import pandas as pd

from src import neighbours as nb


def _gaia(**kw):
    base = dict(source_id=1, ra_2000=10.0, dec_2000=-50.0, phot_g_mean_mag=15.0, bp_rp=np.nan,
                parallax=np.nan, parallax_error=np.nan, teff_gspphot=np.nan, logg_gspphot=np.nan,
                ag_gspphot=np.nan, radius_gspphot=np.nan, radius_flame=np.nan, mass_flame=np.nan)
    base.update(kw)
    return pd.Series(base)


def test_gaia_values_preferred_over_tables():
    p = nb.gaia_parameters(_gaia(bp_rp=1.0, teff_gspphot=5000.0, radius_flame=0.8, mass_flame=0.85))
    assert p["Teff"] == (5000.0, "Gaia GSP-Phot")
    assert p["rad"][1] == "Gaia FLAME" and p["mass"][1] == "Gaia FLAME"


def test_red_dwarf_from_colour_and_parallax():
    # an M dwarf at 64 pc (G = 14.65, BP-RP = 2.88, parallax 16.4 +- 0.4 mas): the target's companion
    p = nb.gaia_parameters(_gaia(phot_g_mean_mag=14.654, bp_rp=2.878, parallax=16.41, parallax_error=0.39))
    assert 3100 < p["Teff"][0] < 3250
    assert 0.25 < p["rad"][0] < 0.45 and p["rad"][1].startswith("G + parallax")
    assert abs(p["mass"][0] - p["rad"][0]) < 0.1  # M ~ R for M dwarfs


def test_bright_red_star_with_tiny_parallax_is_a_giant():
    # TIC 214311906: G = 12.66, BP-RP = 4.68, parallax 0.15 +- 0.09 mas. As a dwarf it would be at 3 pc.
    g = _gaia(phot_g_mean_mag=12.66, bp_rp=4.68, parallax=0.1495, parallax_error=0.0872)
    assert nb.not_a_dwarf(g)
    p = nb.gaia_parameters(g)
    assert "giant" in p["Teff"][1] and 3000 < p["Teff"][0] < 3400
    assert p["rad"][0] > 20 and p["mass"] == (1.0, "evolved star, 1 M_sun assumed")


def test_no_colour_no_estimate():
    assert nb.gaia_parameters(_gaia()) == {}


def test_subluminous_parallax_falls_back_to_dwarf_table():
    p = nb.gaia_parameters(_gaia(phot_g_mean_mag=19.3, bp_rp=1.1, parallax=5.0, parallax_error=1.0))
    assert p["rad"][1].startswith("dwarf table")


def test_match_prefers_magnitude_agreement_over_distance():
    stars = pd.DataFrame({"ID": ["a"], "ra": [10.0], "dec": [-50.0], "Tmag": [13.37]})
    cosd = np.cos(np.radians(-50.0))
    gaia = pd.DataFrame({
        "source_id": [5925373272827253632, 5925373272850782208],
        "ra_2000": [10.0 + 0.8 / 3600 / cosd, 10.0 + 2.2 / 3600 / cosd], "dec_2000": [-50.0, -50.0],
        "phot_g_mean_mag": [15.74, 14.654], "bp_rp": [np.nan, 2.878]})
    ids = nb.match_gaia(stars, gaia)
    assert ids.iloc[0] == 5925373272850782208  # exact 19-digit id survives


def test_fill_parameters_records_sources_and_keeps_tic_values():
    stars = pd.DataFrame({"ID": ["t", "n1", "n2"], "ra": [10.0, 10.0, 10.001], "dec": [-50.0, -50.001, -50.0],
                          "Tmag": [13.0, 14.0, 18.0], "Teff": [3300.0, np.nan, 5000.0],
                          "rad": [0.38, np.nan, 0.8], "mass": [0.37, np.nan, 0.8]})
    gaia = pd.DataFrame([_gaia(source_id=7, ra_2000=10.0, dec_2000=-50.001, phot_g_mean_mag=14.4, bp_rp=1.0)])
    log = nb.fill_parameters(stars, [1, 2], gaia)
    assert len(log) == 1 and log[0]["ID"] == "n1" and log[0]["gaia_dr3"] == 7
    assert np.isfinite(stars.loc[1, ["Teff", "rad", "mass"]].astype(float)).all()
    assert stars.loc[2, "Teff"] == 5000.0  # complete TIC values untouched


def test_load_exclusions(tmp_path):
    f = tmp_path / "pixel_excluded.json"
    assert nb.load_exclusions(f)["tic"] == set()
    f.write_text(json.dumps({"excluded_tic": ["214311906", 214311910],
                             "excluded_gaia_dr3": [5925373268528684288], "reason": {"214311906": "giant"}}))
    e = nb.load_exclusions(f)
    assert e["tic"] == {"214311906", "214311910"} and 5925373268528684288 in e["gaia"]


def test_bc_g_solar():
    assert abs(nb.bc_g(5772.0) - 0.06) < 1e-9
