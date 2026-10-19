"""Neighbour stars for TRICERATOPS: fill missing stellar parameters, drop TIC duplicates, apply exclusions.

TRICERATOPS computes a "planet on this neighbour" (NTP) and "eclipsing binary on this neighbour" (NEB,
NEBx2P) scenario for every catalogued star that could produce the dip. It takes each star's mass, radius
and Teff from the TIC, and when the TIC has none (most faint stars, and many giants) it silently assumes
the Sun: M = R = 1, Teff = 5777 K. For TIC 214311908 that turned a red giant at 17.8" into a Sun-like
star, and "a planet on the giant" became 0.44 of the 0.85 FPP. This module replaces the solar defaults
with values derived from Gaia DR3, in this order of preference:

  Teff: Gaia GSP-Phot Teff -> BP-RP colour through a main-sequence table (Pecaut & Mamajek 2013)
  R:    FLAME radius -> GSP-Phot radius -> Stefan-Boltzmann from G, parallax (S/N >= 3) and a bolometric
        correction (Andrae et al. 2018) -> main-sequence table at that Teff
  M:    FLAME mass -> g R^2 / G from GSP-Phot log g -> 1 M_sun for evolved stars (R > 2 R_sun) ->
        main-sequence table at that Teff

These are rough (no extinction unless GSP-Phot gives A_G; the colour table assumes a dwarf), but much
closer than the Sun for an M dwarf or a giant. Every filled value is recorded with its source.

TIC entries flagged DUPLICATE are a second copy of another star (for TIC 214311908, the co-moving
companion at 9.9" appears twice with Tmag 13.37); they get no flux and no scenarios.

Exclusions: data/followup/<tic>/pixel_excluded.json lists neighbours ruled out as the source of the dip by
a per-pixel PRF fit (src/pixelfit.py). They keep their light in the aperture (they still dilute the
transit) but get no scenarios of their own.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

G_CGS = 6.674e-8
M_SUN_G = 1.989e33
R_SUN_CM = 6.957e10
TEFF_SUN = 5772.0
MBOL_SUN = 4.74
NO_FLUX_TMAG = 40.0  # a star given this magnitude contributes no light to the aperture

# Main-sequence dwarfs: Gaia BP-RP, Teff [K], R [R_sun], M [M_sun], absolute G (Pecaut & Mamajek 2013,
# version 2022; rounded)
DWARFS = pd.DataFrame(
    [(-0.03, 9700, 2.19, 2.18, 1.0), (0.15, 8100, 1.70, 1.86, 1.9), (0.38, 7220, 1.53, 1.61, 2.5),
     (0.59, 6510, 1.36, 1.33, 3.3), (0.76, 5920, 1.10, 1.06, 4.2), (0.85, 5660, 0.95, 0.98, 4.9),
     (0.98, 5280, 0.86, 0.88, 5.6), (1.43, 4440, 0.72, 0.70, 6.9), (1.84, 3850, 0.59, 0.57, 8.2),
     (2.18, 3560, 0.47, 0.44, 9.3), (2.43, 3430, 0.39, 0.37, 10.2), (2.78, 3210, 0.26, 0.23, 11.6),
     (3.24, 3060, 0.20, 0.16, 13.0), (3.80, 2810, 0.14, 0.11, 15.0)],
    columns=["bp_rp", "teff", "rad", "mass", "abs_g"])
# Red giants: Gaia BP-RP -> Teff [K] (approximate, solar metallicity, no reddening)
GIANTS = pd.DataFrame([(0.9, 5100), (1.1, 4750), (1.3, 4500), (1.6, 4150), (2.0, 3850), (2.5, 3600),
                       (3.0, 3450), (4.0, 3300), (5.0, 3150)], columns=["bp_rp", "teff"])

# Andrae et al. (2018) bolometric correction in G, polynomial in (Teff - 5772 K)
_BC_HOT = (6.000e-02, 6.731e-05, -6.647e-08, 2.859e-11, -7.197e-15)  # 4000-8000 K
_BC_COOL = (1.749e00, 1.977e-03, 3.737e-07, -8.966e-11, -4.183e-14)  # 3300-4000 K


def bc_g(teff: float) -> float:
    """Bolometric correction in Gaia G (Andrae et al. 2018), clipped to 3300-8000 K."""
    t = float(np.clip(teff, 3300.0, 8000.0))
    a = _BC_COOL if t < 4000 else _BC_HOT
    x = t - TEFF_SUN
    return float(sum(c * x ** k for k, c in enumerate(a)))


def dwarf(teff: float | None = None, bp_rp: float | None = None) -> dict:
    """Main-sequence Teff, R, M from Teff or (if Teff is unknown) from BP-RP, by linear interpolation."""
    d = DWARFS
    if teff is not None and np.isfinite(teff):
        x, xs = float(np.clip(teff, d.teff.min(), d.teff.max())), d.teff.to_numpy()[::-1]
        return {"Teff": float(teff), "rad": float(np.interp(x, xs, d.rad.to_numpy()[::-1])),
                "mass": float(np.interp(x, xs, d.mass.to_numpy()[::-1]))}
    x = float(np.clip(bp_rp, d.bp_rp.min(), d.bp_rp.max()))
    return {"Teff": float(np.interp(x, d.bp_rp, d.teff)), "rad": float(np.interp(x, d.bp_rp, d.rad)),
            "mass": float(np.interp(x, d.bp_rp, d.mass))}


def _ok(x) -> bool:
    return x is not None and x is not pd.NA and not (isinstance(x, float) and np.isnan(x)) \
        and np.isfinite(float(x))


def gaia_tmag(g_mag: float, bp_rp: float | None) -> float:
    """TESS magnitude from Gaia G and BP-RP (Stassun et al. 2019, eq. 1); G - 0.43 without a colour."""
    if not _ok(bp_rp):
        return float(g_mag) - 0.43
    c = float(np.clip(bp_rp, -0.2, 3.0))
    return float(g_mag) - 0.00522555 * c ** 3 + 0.0891337 * c ** 2 - 0.633923 * c + 0.0324473


def not_a_dwarf(g: pd.Series) -> bool:
    """True when the measured parallax rules out a main-sequence star of this colour: a dwarf this bright
    would be so close that its parallax would be at least twice the measured value + 3 sigma (e.g. a very
    red Tmag-11 star with a 0.15 mas parallax is a distant giant, not an M dwarf at 3 pc)."""
    if not (_ok(g.get("bp_rp")) and _ok(g.get("parallax")) and _ok(g.get("parallax_error"))
            and _ok(g.get("phot_g_mean_mag"))):
        return False
    d = DWARFS
    abs_g = float(np.interp(float(np.clip(g["bp_rp"], d.bp_rp.min(), d.bp_rp.max())), d.bp_rp, d.abs_g))
    plx_dwarf = 1000.0 / 10 ** ((float(g["phot_g_mean_mag"]) - abs_g + 5) / 5)
    return float(g["parallax"]) + 3 * float(g["parallax_error"]) < 0.5 * plx_dwarf


def gaia_parameters(g: pd.Series) -> dict:
    """Best-available Teff, R, M of one Gaia DR3 source (row from `query_gaia`), each with its source."""
    out = {}
    evolved = not_a_dwarf(g)
    if _ok(g.get("teff_gspphot")):
        out["Teff"] = (float(g["teff_gspphot"]), "Gaia GSP-Phot")
    elif _ok(g.get("bp_rp")) and evolved:
        x = float(np.clip(g["bp_rp"], GIANTS.bp_rp.min(), GIANTS.bp_rp.max()))
        out["Teff"] = (float(np.interp(x, GIANTS.bp_rp, GIANTS.teff)), "BP-RP, giant table (parallax rules out a dwarf)")
    elif _ok(g.get("bp_rp")):
        out["Teff"] = (dwarf(bp_rp=float(g["bp_rp"]))["Teff"], "BP-RP, dwarf table")
    else:
        return out  # without a colour nothing else can be estimated
    teff = out["Teff"][0]
    if _ok(g.get("radius_flame")):
        out["rad"] = (float(g["radius_flame"]), "Gaia FLAME")
    elif _ok(g.get("radius_gspphot")):
        out["rad"] = (float(g["radius_gspphot"]), "Gaia GSP-Phot")
    elif _ok(g.get("parallax")) and _ok(g.get("parallax_error")) and _ok(g.get("phot_g_mean_mag")) and (
            evolved or (g["parallax"] > 0 and g["parallax"] / g["parallax_error"] >= 3)):
        # for a star that cannot be a dwarf, a low-S/N parallax still gives the size to a factor ~2;
        # it is floored at 1 sigma so a negative or tiny parallax does not make the star absurdly far
        plx = max(float(g["parallax"]), float(g["parallax_error"]))
        ag = float(g["ag_gspphot"]) if _ok(g.get("ag_gspphot")) else 0.0
        m_g = float(g["phot_g_mean_mag"]) + 5 * np.log10(plx / 100.0) - ag
        lum = 10 ** (-0.4 * (m_g + bc_g(teff) - MBOL_SUN))
        out["rad"] = (float(np.sqrt(lum) / (teff / TEFF_SUN) ** 2),
                      "G + parallax + BC" + ("" if ag else " (no extinction)")
                      + (", low-S/N parallax" if g["parallax"] / g["parallax_error"] < 3 else ""))
        if not evolved and out["rad"][0] < 0.5 * dwarf(teff=teff)["rad"]:
            # far below the main sequence: a white dwarf or a bad parallax; TRICERATOPS models neither
            out["rad"] = (dwarf(teff=teff)["rad"], "dwarf table (parallax implies a subluminous star)")
    else:
        out["rad"] = (dwarf(teff=teff)["rad"], "dwarf table")
    rad = out["rad"][0]
    if _ok(g.get("mass_flame")):
        out["mass"] = (float(g["mass_flame"]), "Gaia FLAME")
    elif _ok(g.get("logg_gspphot")):
        m = 10 ** float(g["logg_gspphot"]) * (rad * R_SUN_CM) ** 2 / G_CGS / M_SUN_G
        out["mass"] = (float(np.clip(m, 0.08, 10.0)), "GSP-Phot log g and radius")
    elif rad > 2.0 or evolved:
        out["mass"] = (1.0, "evolved star, 1 M_sun assumed")
    elif out["rad"][1].startswith("G + parallax"):
        d = DWARFS.sort_values("rad")
        out["mass"] = (float(np.interp(rad, d.rad, d.mass)), "dwarf mass-radius relation")
    else:
        out["mass"] = (dwarf(teff=teff)["mass"], "dwarf table")
    return out


def match_gaia(stars: pd.DataFrame, gaia: pd.DataFrame, radius_arcsec: float = 3.0,
               dmag_max: float = 1.5) -> pd.Series:
    """Gaia source_id for each star (index of `stars`), matched at epoch 2000, or <NA>.

    TIC positions are epoch J2000; `gaia` must carry ra_2000/dec_2000 (Gaia propagated to 2000). Old
    2MASS-based TIC entries can sit ~2" from their Gaia counterpart, so the radius is generous and the
    match minimises (distance / 1")^2 + (Tmag difference / 0.5 mag)^2, with Tmag predicted from G and
    BP-RP. (For TIC 214311908, the companion's 2MASS entry is 0.8" from a fainter Gaia source but 2.2"
    from the companion itself, whose predicted Tmag agrees to 0.1 mag.)"""
    out = pd.Series(pd.NA, index=stars.index, dtype="Int64")  # 19-digit ids do not survive float
    if gaia.empty:
        return out
    cosd = np.cos(np.radians(stars["dec"].to_numpy()))
    for k, (i, s) in enumerate(stars.iterrows()):
        d = np.hypot((gaia["ra_2000"] - s["ra"]) * cosd[k], gaia["dec_2000"] - s["dec"]) * 3600
        tmag = np.array([gaia_tmag(gm, c) if _ok(gm) else np.nan
                         for gm, c in zip(gaia["phot_g_mean_mag"], gaia.get("bp_rp", [np.nan] * len(gaia)))])
        dm = np.abs(tmag - float(s["Tmag"]))
        ok = (d < radius_arcsec) & (dm < dmag_max)
        if ok.any():
            score = np.where(ok, (d / 1.0) ** 2 + (dm / 0.5) ** 2, np.inf)
            out[i] = int(gaia["source_id"].iloc[int(np.argmin(score))])
    return out


def fill_parameters(stars: pd.DataFrame, rows: list, gaia: pd.DataFrame) -> list[dict]:
    """Replace missing Teff/rad/mass of stars.loc[rows] in place with Gaia-derived values.

    Returns one record per star touched: TIC ID, Gaia source, the filled values and their sources, or
    the reason nothing could be filled (TRICERATOPS then falls back to the Sun for that star)."""
    ids = match_gaia(stars.loc[rows], gaia)
    log = []
    for i in rows:
        miss = [p for p in ("Teff", "rad", "mass") if not _ok(stars.at[i, p])]
        if not miss:
            continue
        rec = {"ID": str(stars.at[i, "ID"]), "Tmag": float(stars.at[i, "Tmag"]),
               "gaia_dr3": None if pd.isna(ids[i]) else int(ids[i])}
        if pd.isna(ids[i]):
            rec["unfilled"] = miss
            rec["why"] = "no Gaia DR3 match within 1.5 arcsec"
        else:
            par = gaia_parameters(gaia.set_index("source_id").loc[int(ids[i])])
            # a TIC radius (giants often have R and Teff but no mass) is kept and used for the mass
            if "rad" not in miss and "mass" in miss and "Teff" in par:
                rad = float(stars.at[i, "rad"])
                par["mass"] = (1.0, "evolved star, 1 M_sun assumed") if rad > 2.0 else \
                    (dwarf(teff=float(stars.at[i, "Teff"]) if _ok(stars.at[i, "Teff"]) else par["Teff"][0])["mass"],
                     "dwarf table")
            for p in miss:
                if p in par:
                    stars.at[i, p] = par[p][0]
                    rec[p] = [round(par[p][0], 3), par[p][1]]
            left = [p for p in miss if p not in par]
            if left:
                rec["unfilled"] = left
                rec["why"] = "Gaia source has no BP-RP colour"
        log.append(rec)
    return log


def query_gaia(ra: float, dec: float, radius_arcsec: float) -> pd.DataFrame:
    """Gaia DR3 sources in a cone with the columns `gaia_parameters` uses, positions also at J2000."""
    from astroquery.vizier import Vizier
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    c = SkyCoord(ra, dec, unit="deg")
    r = radius_arcsec * u.arcsec
    main = Vizier(columns=["Source", "RAJ2000", "DEJ2000", "Gmag", "BP-RP", "Plx", "e_Plx", "Teff", "logg",
                           "AG"], row_limit=-1).query_region(c, radius=r, catalog="I/355/gaiadr3")
    if not len(main):
        return pd.DataFrame(columns=["source_id", "ra_2000", "dec_2000", "phot_g_mean_mag"])
    df = main[0].to_pandas().rename(columns={
        "Source": "source_id", "RAJ2000": "ra_2000", "DEJ2000": "dec_2000", "Gmag": "phot_g_mean_mag",
        "BP-RP": "bp_rp", "Plx": "parallax", "e_Plx": "parallax_error", "Teff": "teff_gspphot",
        "logg": "logg_gspphot", "AG": "ag_gspphot"})
    par = Vizier(columns=["Source", "Rad", "Rad-Flame", "Mass-Flame"], row_limit=-1).query_region(
        c, radius=r, catalog="I/355/paramp")
    if len(par):
        p = par[0].to_pandas().rename(columns={"Source": "source_id", "Rad": "radius_gspphot",
                                               "Rad-Flame": "radius_flame", "Mass-Flame": "mass_flame"})
        df = df.merge(p, on="source_id", how="left")
    return df


def load_exclusions(path: Path) -> dict:
    """pixel_excluded.json -> {"tic": set of TIC IDs (str), "gaia": set of Gaia DR3 ids, "reason": {...}}."""
    if not path.exists():
        return {"tic": set(), "gaia": set(), "reason": {}}
    d = json.loads(path.read_text())
    return {"tic": {str(x) for x in d.get("excluded_tic", [])},
            "gaia": {int(x) for x in d.get("excluded_gaia_dr3", [])},
            "reason": d.get("reason", {})}
