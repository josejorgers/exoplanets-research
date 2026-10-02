"""Step 6 — cross-match the vetted signals against known objects: planets, TOIs, CTOIs, EBs, neighbours.

Usage:
    uv run python -m src.known                         # sector from .env, catalogs downloaded once per day
    uv run python -m src.known --offline               # use the newest catalogs already in data/known/
    uv run python -m src.known --refresh               # force a new download today
    uv run python -m src.known --gaia-scope strong     # Gaia EB cones around every SDE>=9 signal, not only candidates
    uv run python -m src.known --second-pass           # the second-pass signals (data/search/pass2/)

Input : data/search/sector_{N}_vetted.parquet (from src.vet)
        data/targets/sector_{N}.parquet         (ra, dec, tmag, gaia_dr2_id)
        catalogs, stored raw in data/known/YYYY-MM-DD/:
          toi.csv, ctoi.csv (ExoFOP), pscomppars.csv (NASA Exoplanet Archive),
          tess_ebs.parquet (Prsa+ 2022, VizieR J/ApJS/258/16, sectors 1-26),
          gaia_veb_sector_{N}.parquet + _queried.json (Gaia DR3 EBs, VizieR I/358/veb, cone queries)
Output: data/search/sector_{N}_known.parquet   all vetted columns + known_class and match details.
        Sorted: passes_all first, then known_class (new first), then snr.
        data/search/sector_{N}_known_meta.json catalog folder, sizes, tolerances.

known_class (first rule that applies): known_planet, toi, ctoi, known_eb (same star, same ephemeris),
neighbour (a catalogued object within --radius has this ephemeris), sibling (another star of this
sector has the same signal and can explain it), new_in_known_system (the star has catalogued
objects, none at this period), new.

Walkthrough of every step, with worked examples: notebooks/pipeline-6-known-walkthrough.ipynb
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
import warnings
from datetime import date
from pathlib import Path

import astropy.units as u
import dotenv
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier

log = logging.getLogger("known")

dotenv.load_dotenv()
DATA_DIR = Path(
    os.environ.get("DATA_DIR", Path(__file__).resolve().parents[1] / "data")
)
SEARCH_DIR = DATA_DIR / "search"
KNOWN_ROOT = DATA_DIR / "known"

BTJD0 = 2457000.0  # BTJD = BJD - 2457000
RATIOS = (1 / 3, 1 / 2, 1, 2, 3)  # P_ours / P_catalog aliases accepted
PERIOD_TOL = 0.002  # relative; real matches agree to < 0.05 %, chance match ~0.5 % per object
NEIGHBOUR_RADIUS_ARCSEC = 120.0  # ~6 TESS pixels (21")
SELF_SEP_ARCSEC = 3.0  # Gaia EB counts as the target itself (fallback when the Gaia IDs differ)
CLASS_ORDER = ["new", "new_in_known_system", "sibling", "neighbour", "known_eb", "ctoi", "toi", "known_planet"]

URLS = {
    "toi.csv": "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv",
    "ctoi.csv": "https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv",
    "pscomppars.csv": (
        "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query="
        "select+pl_name,hostname,tic_id,ra,dec,sy_tmag,pl_orbper,pl_orbpererr1,pl_tranmid,pl_tranmiderr1,"
        "pl_trandep,pl_rade,tran_flag+from+pscomppars"
    ),
}
KNOWN_COLS = [
    "catalog", "obj_id", "name", "tic_id", "ra", "dec", "tmag", "period", "period_err",
    "epoch", "epoch_err", "depth_ppm", "disposition",
]


# ------------------------------------------------------------------------------------------------
# catalogs
# ------------------------------------------------------------------------------------------------


def catalog_dir(offline: bool = False) -> Path:
    """Today's folder, or (offline) the newest existing one."""
    if offline:
        days = sorted(p for p in KNOWN_ROOT.glob("????-??-??") if (p / "toi.csv").exists())
        if not days:
            raise FileNotFoundError(f"--offline but no catalogs in {KNOWN_ROOT}")
        return days[-1]
    d = KNOWN_ROOT / date.today().isoformat()
    d.mkdir(parents=True, exist_ok=True)
    return d


def fetch_catalogs(day_dir: Path, refresh: bool = False, offline: bool = False) -> dict[str, pd.DataFrame]:
    """Download (once per day) and read the four static catalogs."""
    for name, url in URLS.items():
        path = day_dir / name
        if offline or (path.exists() and not refresh):
            continue
        log.info("downloading %s", name)
        r = requests.get(url, timeout=300)
        r.raise_for_status()
        path.write_bytes(r.content)
    eb_path = day_dir / "tess_ebs.parquet"
    if not offline and (refresh or not eb_path.exists()):
        log.info("downloading TESS EB catalog (VizieR J/ApJS/258/16)")
        t = Vizier(
            columns=["TIC", "m_TIC", "RAJ2000", "DEJ2000", "Tmag", "BJD0", "e_BJD0", "Per", "e_Per", "Morph"],
            row_limit=-1,
        ).get_catalogs("J/ApJS/258/16/tess-ebs")[0]
        t.to_pandas().to_parquet(eb_path, index=False)
    return {
        "toi": pd.read_csv(day_dir / "toi.csv"),
        "ctoi": pd.read_csv(day_dir / "ctoi.csv"),
        "planet": pd.read_csv(day_dir / "pscomppars.csv"),
        "tess_eb": pd.read_parquet(eb_path),
    }


def _btjd(bjd):
    """Full BJD -> BTJD; anything outside the TESS era (typos in community tables) -> NaN."""
    bjd = pd.to_numeric(bjd, errors="coerce")
    return (bjd - BTJD0).where((bjd > 2458000) & (bjd < 2470000))


def _norm_toi(t: pd.DataFrame) -> pd.DataFrame:
    c = SkyCoord(t["RA"].values, t["Dec"].values, unit=(u.hourangle, u.deg))  # sexagesimal strings
    return pd.DataFrame({
        "catalog": "toi", "obj_id": t["TOI"].map(lambda x: f"TOI-{x:.2f}"),
        "name": t["Planet Name"].fillna(""), "tic_id": t["TIC ID"], "ra": c.ra.deg, "dec": c.dec.deg,
        "tmag": t["TESS Mag"], "period": t["Period (days)"], "period_err": t["Period (days) err"],
        "epoch": _btjd(t["Epoch (BJD)"]), "epoch_err": t["Epoch (BJD) err"], "depth_ppm": t["Depth (ppm)"],
        "disposition": t["TFOPWG Disposition"].fillna(t["TESS Disposition"]),
    })


def _norm_ctoi(t: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({
        "catalog": "ctoi", "obj_id": t["CTOI"].map(lambda x: f"CTOI-{x:.2f}"),
        "name": t["Promoted to TOI"].map(lambda x: f"TOI-{x:.2f}" if pd.notna(x) else ""),
        "tic_id": t["TIC ID"], "ra": t["RA"], "dec": t["Dec"], "tmag": t["TESS Mag"],
        "period": t["Period (days)"], "period_err": t["Period (days) Error"],
        "epoch": _btjd(t["Transit Epoch (BJD)"]), "epoch_err": t["Transit Epoch (BJD) err"],
        "depth_ppm": t["Depth ppm"], "disposition": t["TFOPWG Disposition"].fillna(t["User Disposition"]),
    })


def _norm_planet(t: pd.DataFrame) -> pd.DataFrame:
    t = t[t["tic_id"].notna()]
    return pd.DataFrame({
        "catalog": "planet", "obj_id": t["pl_name"], "name": t["pl_name"],
        "tic_id": t["tic_id"].str.replace("TIC ", "").astype("int64"), "ra": t["ra"], "dec": t["dec"],
        "tmag": t["sy_tmag"], "period": t["pl_orbper"], "period_err": t["pl_orbpererr1"],
        "epoch": _btjd(t["pl_tranmid"]), "epoch_err": t["pl_tranmiderr1"],
        "depth_ppm": t["pl_trandep"] * 1e4,  # the archive gives depth in %
        "disposition": "CP",
    })


def _norm_tess_eb(t: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({
        "catalog": "tess_eb",
        "obj_id": "TESS-EB " + t["TIC"].astype(str).str.lstrip("0") + "." + t["m_TIC"].astype(str),
        "name": "", "tic_id": t["TIC"].astype("int64"), "ra": t["RAJ2000"], "dec": t["DEJ2000"],
        "tmag": t["Tmag"], "period": t["Per"], "period_err": t["e_Per"],
        "epoch": t["BJD0"], "epoch_err": t["e_BJD0"],  # already BTJD in this catalog
        "depth_ppm": np.nan, "disposition": "EB",
    })


def load_known(raw: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """All static catalogs in one schema (KNOWN_COLS)."""
    known = pd.concat(
        [_norm_toi(raw["toi"]), _norm_ctoi(raw["ctoi"]), _norm_planet(raw["planet"]), _norm_tess_eb(raw["tess_eb"])],
        ignore_index=True,
    )
    known["tic_id"] = known["tic_id"].astype("Int64")
    for c in ["ra", "dec", "tmag", "period", "period_err", "epoch", "epoch_err", "depth_ppm"]:
        known[c] = pd.to_numeric(known[c], errors="coerce")
    return known[KNOWN_COLS]


def gaia_ebs(sig: pd.DataFrame, day_dir: Path, sector: int, offline: bool) -> pd.DataFrame:
    """Gaia DR3 EBs within the neighbour radius of the signals in `sig` (VizieR cone queries, cached).

    The cache records WHICH signals were queried (a cone with no EB has no rows), so a changed
    candidate list only queries the new ones. Online: today's folder, missing cones are queried.
    Offline: the caches of all dated folders are merged; cones still missing are reported.
    Returned in the KNOWN_COLS schema plus tic_query (the signal it was found around) and self.
    """
    cols = ["tic_query", "Source", "RA_ICRS", "DE_ICRS", "Freq"]
    dirs = sorted(KNOWN_ROOT.glob("????-??-??"), reverse=True) if offline else [day_dir]
    parts, queried = [pd.DataFrame(columns=cols)], set()
    for d in dirs:
        rows_p, q_p = d / f"gaia_veb_sector_{sector}.parquet", d / f"gaia_veb_sector_{sector}_queried.json"
        if rows_p.exists() and q_p.exists():
            new_q = set(json.loads(q_p.read_text())) - queried
            part = pd.read_parquet(rows_p)
            parts.append(part[part["tic_query"].isin(new_q)])
            queried |= new_q
    parts = [p for p in parts if len(p)] or parts[:1]
    g = pd.concat(parts, ignore_index=True)
    missing = sig[~sig.index.isin(queried)]
    if len(missing) and offline:
        log.warning("offline: no cached Gaia EB query for %d of %d signals; they get no Gaia check",
                    len(missing), len(sig))
    elif len(missing):
        log.info("querying Gaia DR3 EBs (VizieR I/358/veb) around %d signals", len(missing))
        res = Vizier(columns=["Source", "RA_ICRS", "DE_ICRS", "Freq"], row_limit=-1, timeout=600).query_region(
            SkyCoord(missing["ra"].values, missing["dec"].values, unit="deg"),
            radius=NEIGHBOUR_RADIUS_ARCSEC * u.arcsec,
            catalog="I/358/veb",
        )
        new = res[0].to_pandas() if len(res) else pd.DataFrame(columns=["_q"] + cols[1:])
        new["tic_query"] = missing.index.values[new["_q"].astype(int) - 1]  # _q is 1-based
        g = pd.concat([x for x in (g, new[cols]) if len(x)] or [g], ignore_index=True)
        queried |= set(int(t) for t in missing.index)
        g.to_parquet(day_dir / f"gaia_veb_sector_{sector}.parquet", index=False)
        (day_dir / f"gaia_veb_sector_{sector}_queried.json").write_text(json.dumps(sorted(int(t) for t in queried)))
    g = g[g["tic_query"].isin(sig.index)].reset_index(drop=True)
    if len(g) == 0:
        return pd.DataFrame(columns=KNOWN_COLS + ["tic_query", "self"])
    out = pd.DataFrame({
        "catalog": "gaia_eb", "obj_id": "Gaia DR3 " + g["Source"].astype(str), "name": "",
        "tic_id": pd.array([pd.NA] * len(g), dtype="Int64"), "ra": g["RA_ICRS"], "dec": g["DE_ICRS"],
        "tmag": np.nan, "period": 1 / g["Freq"], "period_err": np.nan, "epoch": np.nan, "epoch_err": np.nan,
        "depth_ppm": np.nan, "disposition": "EB", "tic_query": g["tic_query"].astype("int64"),
    })
    host = sig.reindex(out["tic_query"])
    sep = SkyCoord(out["ra"].values, out["dec"].values, unit="deg").separation(
        SkyCoord(host["ra"].values, host["dec"].values, unit="deg")
    ).arcsec
    # Gaia DR3 source IDs are (almost always) the DR2 IDs the TIC stores; positions differ by epoch (2016 vs 2000)
    same_id = g["Source"].astype("int64").values == host["gaia_dr2_id"].fillna(-1).astype("int64").values
    out["self"] = same_id | (sep < SELF_SEP_ARCSEC)
    return out.reset_index(drop=True)


# ------------------------------------------------------------------------------------------------
# matching
# ------------------------------------------------------------------------------------------------


def period_ratio(p_ours: float, p_cat: float, tol: float = PERIOD_TOL):
    """(r, relative difference) for the alias r with p_ours ~ r * p_cat, or None."""
    if not (np.isfinite(p_ours) and np.isfinite(p_cat)) or p_cat <= 0:
        return None
    r, d = min(((r, abs(p_ours - r * p_cat) / p_ours) for r in RATIOS), key=lambda x: x[1])
    return (r, d) if d < tol else None


def ephemeris_check(T0: float, dur_h: float, k, r: float) -> tuple[float, float, str]:
    """Our epoch vs catalog row k's ephemeris, for r = P_ours / P_cat: (dt [d], sigma_t [d], ok|bad|unknown).

    Propagates with the CATALOG period (fitted over years, far more precise than one sector's):
    propagating with ours turns a 0.002-d period difference into hours over hundreds of orbits.
    """
    if not (np.isfinite(k.epoch) and np.isfinite(k.period)):
        return np.nan, np.nan, "unknown"
    m = min(r, 1) * k.period  # spacing of the predicted events
    n = round((T0 - k.epoch) / m)
    e_err = k.epoch_err if np.isfinite(k.epoch_err) else np.inf  # missing error -> prediction unusable
    p_err = k.period_err if np.isfinite(k.period_err) else np.inf
    sigma_t = float(np.hypot(e_err, n * min(r, 1) * p_err))
    dt = T0 - (k.epoch + n * m)
    if k.catalog in ("tess_eb", "gaia_eb"):  # our dip may be the secondary eclipse
        dt2 = (T0 - k.epoch) % m - m / 2  # offset from phase 0.5
        dt = dt if abs(dt) <= abs(dt2) else dt2
    if not np.isfinite(sigma_t) or sigma_t > m / 4:
        return dt, sigma_t, "unknown"
    return dt, sigma_t, "ok" if abs(dt) < max(3 * sigma_t, dur_h / 24) else "bad"


def match_row(s, k) -> dict | None:
    """Period + ephemeris match of signal s against catalog row k; None if no match."""
    pr = period_ratio(s.period, k.period)
    if pr is None:
        return None
    dt, st, status = ephemeris_check(s.T0, s.duration_h, k, pr[0])
    if status == "bad":
        return None
    return {"ratio": pr[0], "rel_diff": pr[1], "dt_min": dt * 1440, "ephemeris": status}


def depth_needed_ppm(depth_obs_ppm: float, tmag_target: float, tmag_nb: float) -> float:
    """Minimum eclipse depth a neighbour needs to produce our observed depth.

    Our light curves are deblended (SPOC CROWDSAP; QLP shown in the walkthrough §6), which turns a
    neighbour's eclipse D_n into D_obs ~ f D_n F_n / F_t, f <= 1 -> D_n >= D_obs F_t / F_n.
    """
    return depth_obs_ppm * 10 ** (-0.4 * (tmag_target - tmag_nb))


def likely_source(depth_obs_ppm, tmag_target, tmag_nb, depth_nb_ppm) -> tuple[float, str]:
    need = depth_needed_ppm(depth_obs_ppm, tmag_target, tmag_nb)
    if not np.isfinite(need):
        return need, "unclear"
    if need > 1e6:
        return need, "target"  # the neighbour is too faint to make our dip at all
    if not np.isfinite(depth_nb_ppm):
        return need, "unclear"
    ratio = depth_nb_ppm / need
    return need, "neighbour" if ratio >= 0.7 else "target" if ratio < 0.3 else "unclear"


def find_siblings(sig: pd.DataFrame, radius_arcsec: float) -> pd.DataFrame:
    """Pairs of strong signals (SDE and transit count pass) in this sector with the same ephemeris."""
    strong = sig[sig["pass_sde"] & sig["pass_n_transits"] & sig["ra"].notna()]
    cols = ["tic_a", "tic_b", "sep_arcsec", "dt_min"]
    if len(strong) < 2:
        return pd.DataFrame(columns=cols)
    sc = SkyCoord(strong["ra"].values, strong["dec"].values, unit="deg")
    i, j, sep, _ = sc.search_around_sky(sc, radius_arcsec * u.arcsec)
    rows = []
    for a, b, d in zip(i, j, sep):
        if a >= b:
            continue
        sa, sb = strong.iloc[a], strong.iloc[b]
        if period_ratio(sa.period, sb.period) is None:
            continue
        m = min(sa.period, sb.period)
        dt = (sa.T0 - sb.T0 + m / 2) % m - m / 2
        if abs(dt) * 24 > max(sa.duration_h, sb.duration_h):  # same sector: epochs must line up
            continue
        rows.append((sa.name, sb.name, d.arcsec, dt * 1440))
    return pd.DataFrame(rows, columns=cols)


def classify_all(sig: pd.DataFrame, known: pd.DataFrame, gaia: pd.DataFrame, siblings: pd.DataFrame,
                 radius_arcsec: float) -> pd.DataFrame:
    """One row of match results per signal (index tic_id)."""
    # all (signal, catalog row) pairs within the radius, in one k-d tree search
    has_pos = sig["ra"].notna()
    sc = SkyCoord(sig.loc[has_pos, "ra"].values, sig.loc[has_pos, "dec"].values, unit="deg")
    kc = SkyCoord(known["ra"].values, known["dec"].values, unit="deg")
    ik, isig, sep, _ = sc.search_around_sky(kc, radius_arcsec * u.arcsec)  # ik -> known, isig -> sc
    near = pd.DataFrame({"tic": sig.index[has_pos][isig], "k": ik, "sep": sep.arcsec})
    near = near[known["tic_id"].fillna(-1).values[near["k"]] != near["tic"].values]  # other stars only
    near_by_tic = {t: g for t, g in near.groupby("tic")}
    same_by_tic = {t: g for t, g in known[known["tic_id"].isin(sig.index)].groupby("tic_id")}
    gaia_by_tic = {t: g for t, g in gaia.groupby("tic_query")} if len(gaia) else {}
    flag_of = {"planet": "known_planet", "toi": "known_toi", "ctoi": "known_ctoi", "tess_eb": "known_eb"}

    results = {}
    for tic, s in sig.iterrows():
        out = {
            "known_planet": False, "known_toi": False, "known_ctoi": False, "known_eb": False,
            "neighbour_match": False, "sibling_match": False, "system_known": False, "known_ids": "",
            "nb_id": "", "nb_sep_arcsec": np.nan, "nb_depth_needed_ppm": np.nan, "likely_source": "",
        }
        ids = []
        # --- same star ---
        same = same_by_tic.get(tic)
        if same is not None:
            out["system_known"] = True
            for k in same.itertuples():
                if match_row(s, k):
                    out[flag_of[k.catalog]] = True
                    ids.append(f"{k.obj_id} [{k.disposition}]")
        g = gaia_by_tic.get(tic)
        if g is not None:
            for k in g[g["self"]].itertuples():
                if period_ratio(s.period, k.period):
                    out["known_eb"] = True
                    ids.append(f"{k.obj_id} [EB]")
        # --- neighbours: catalogs and Gaia; keep the closest match ---
        best = []
        nb = near_by_tic.get(tic)
        if nb is not None:
            for kidx, d in zip(nb["k"], nb["sep"]):
                k = known.iloc[kidx]
                if match_row(s, k):
                    need, src = likely_source(s.depth_ppm, s.tmag, k.tmag, k.depth_ppm)
                    best.append((d, f"{k.obj_id} on TIC {k.tic_id} [{k.disposition}]", need, src))
        if g is not None:
            gsep = SkyCoord(g["ra"].values, g["dec"].values, unit="deg").separation(
                SkyCoord(s.ra, s.dec, unit="deg")).arcsec
            for k, d in zip(g.itertuples(), gsep):
                if not k.self and period_ratio(s.period, k.period):
                    best.append((d, f"{k.obj_id} [EB]", np.nan, "unclear"))
        if best:
            out["neighbour_match"] = True
            out["nb_sep_arcsec"], out["nb_id"], out["nb_depth_needed_ppm"], out["likely_source"] = min(
                best, key=lambda x: x[0]
            )
        # --- siblings in this sector: the other star must be able to explain our dip ---
        if len(siblings):
            pairs = siblings[(siblings["tic_a"] == tic) | (siblings["tic_b"] == tic)]
            for r in pairs.itertuples():
                other = sig.loc[r.tic_b if r.tic_a == tic else r.tic_a]
                need = depth_needed_ppm(s.depth_ppm, s.tmag, other.tmag)
                if need <= 1.43 * other.depth_ppm:  # same 0.7 margin as likely_source
                    out["sibling_match"] = True
                    if not best:
                        out["nb_id"], out["nb_sep_arcsec"] = f"TIC {other.name} (this sector)", r.sep_arcsec
                        out["nb_depth_needed_ppm"], out["likely_source"] = need, "neighbour"
        out["known_ids"] = "; ".join(ids)
        out["known_class"] = (
            "known_planet" if out["known_planet"] else "toi" if out["known_toi"]
            else "ctoi" if out["known_ctoi"] else "known_eb" if out["known_eb"]
            else "neighbour" if out["neighbour_match"] else "sibling" if out["sibling_match"]
            else "new_in_known_system" if out["system_known"] else "new"
        )
        results[tic] = out
    return pd.DataFrame.from_dict(results, orient="index")


# ------------------------------------------------------------------------------------------------
# sector
# ------------------------------------------------------------------------------------------------


def known_sector(
    sector: int,
    offline: bool = False,
    refresh: bool = False,
    gaia_scope: str = "candidates",
    radius_arcsec: float = NEIGHBOUR_RADIUS_ARCSEC,
    second_pass: bool = False,
) -> pd.DataFrame:
    base = SEARCH_DIR / "pass2" / f"sector_{sector}_pass2" if second_pass else SEARCH_DIR / f"sector_{sector}"
    vetted = pd.read_parquet(f"{base}_vetted.parquet")
    targets = pd.read_parquet(
        DATA_DIR / "targets" / f"sector_{sector}.parquet", columns=["tic_id", "ra", "dec", "tmag", "gaia_dr2_id"]
    )
    sig = vetted.merge(targets, on="tic_id", how="left").set_index("tic_id")
    sig = sig[sig["status"] != "error"]

    day_dir = catalog_dir(offline)
    raw = fetch_catalogs(day_dir, refresh, offline)
    known = load_known(raw)
    log.info("catalogs from %s: %s", day_dir, known["catalog"].value_counts().to_dict())

    scope_rows = {
        "candidates": sig["passes_all"],
        "strong": sig["pass_sde"] & sig["pass_n_transits"],
        "none": pd.Series(False, index=sig.index),
    }[gaia_scope]
    gaia = gaia_ebs(sig[scope_rows & sig["ra"].notna()], day_dir, sector, offline)
    log.info("Gaia DR3 EBs within %.0f\" of %d signals: %d (%d are the target itself)",
             radius_arcsec, int(scope_rows.sum()), len(gaia), int(gaia["self"].sum()) if len(gaia) else 0)
    siblings = find_siblings(sig, radius_arcsec)
    log.info("sibling pairs in sector %d: %d", sector, len(siblings))

    res = classify_all(sig, known, gaia, siblings, radius_arcsec)
    out = sig.drop(columns=["ra", "dec", "tmag", "gaia_dr2_id"]).join(res).reset_index()
    out = out.join(sig[["ra", "dec", "tmag"]], on="tic_id")
    out["_order"] = out["known_class"].map(CLASS_ORDER.index)
    out = out.sort_values(["passes_all", "_order", "snr"], ascending=[False, True, False]).drop(columns="_order")
    out = out.reset_index(drop=True)

    out_path = Path(f"{base}_known.parquet")
    out.to_parquet(out_path, index=False)
    meta = {
        "catalog_dir": str(day_dir),
        "catalog_rows": {k: int(v) for k, v in known["catalog"].value_counts().items()},
        "gaia_scope": gaia_scope,
        "gaia_rows": int(len(gaia)),
        "radius_arcsec": radius_arcsec,
        "period_tol": PERIOD_TOL,
        "period_ratios": list(RATIOS),
        "written": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    out_path.with_name(out_path.stem + "_meta.json").write_text(json.dumps(meta, indent=2))
    _log_summary(out, out_path)
    return out


def _log_summary(out: pd.DataFrame, out_path: Path) -> None:
    cand = out[out["passes_all"]]
    log.info("wrote %s: %d signals, %d pass vetting", out_path, len(out), len(cand))
    log.info("classes among candidates: %s", cand["known_class"].value_counts().to_dict())
    log.info("classes among all signals: %s", out["known_class"].value_counts().to_dict())
    for r in cand[cand["known_class"].isin(["new", "new_in_known_system"])].itertuples():
        log.info(
            "  %-19s TIC %10d  P %8.4f d  depth %7.0f ppm  Rp %5.1f Re  SNR %5.1f  flags: %s",
            r.known_class, r.tic_id, r.period, r.depth_ppm, r.rp_earth, r.snr, r.flags or "-",
        )
    # recall check on vetting: catalogued planets/TOIs whose signal we found but vetting rejected
    missed = out[~out["passes_all"] & out["known_class"].isin(["known_planet", "toi"])]
    log.info("known planets/TOIs recovered by the search but rejected by vetting: %d", len(missed))
    for r in missed.itertuples():
        log.info("  TIC %10d  %-40s failed: %s", r.tic_id, r.known_ids[:40], r.failed)


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--sector", type=int, default=int(os.environ.get("SECTOR", 0)) or None
    )
    p.add_argument("--offline", action="store_true", help="use the newest catalogs already in data/known/")
    p.add_argument("--refresh", action="store_true", help="download the catalogs again even if today's exist")
    p.add_argument(
        "--gaia-scope", choices=["candidates", "strong", "none"], default="candidates",
        help="Gaia EB cones around passes_all signals (default), every SDE>=9 signal with >=2 transits, or none",
    )
    p.add_argument("--radius", type=float, default=NEIGHBOUR_RADIUS_ARCSEC, help="neighbour radius, arcsec")
    p.add_argument("--second-pass", action="store_true",
                   help="cross-match the second-pass signals (pass2/sector_N_pass2_vetted.parquet)")
    a = p.parse_args()
    if a.sector is None:
        p.error("no --sector given and SECTOR not set in .env")
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    warnings.filterwarnings("ignore", module="astroquery")
    out = known_sector(a.sector, a.offline, a.refresh, a.gaia_scope, a.radius, a.second_pass)
    cols = ["tic_id", "known_class", "period", "depth_ppm", "rp_earth", "snr", "flags", "known_ids", "nb_id",
            "likely_source"]
    print(out[out["passes_all"]][cols].to_string())


if __name__ == "__main__":
    main()
