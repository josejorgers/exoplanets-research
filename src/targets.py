"""Step 1 — build the target list of nearby M dwarfs observed in one TESS sector.

Usage:
    uv run python -m src.targets --sector 104
    uv run python -m src.targets --sector 104 --tmag-max 13.5 --force

Output: data/targets/sector_{N}.parquet (one row per TIC star) and
        data/targets/tic_mdwarfs_allsky.parquet (cached all-sky TIC query).
"""

from __future__ import annotations
import dotenv

dotenv.load_dotenv()
import os

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from astroquery.mast import Catalogs
from tess_stars2px import tess_stars2px_function_entry

log = logging.getLogger("targets")

DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "targets"

# Columns kept from the TIC. Full list: Stassun et al. 2019, arXiv:1905.10694, Table 1,
# or `Catalogs.query_criteria(...).colnames`.
TIC_COLUMNS = [
    "ID",
    "ra",
    "dec",
    "Tmag",
    "Teff",
    "rad",
    "mass",
    "d",
    "plx",
    "contratio",
    "lumclass",
    "objType",
    "GAIA",
    "pmRA",
    "pmDEC",
]


def query_tic_mdwarfs(
    teff_max: float = 4000.0,
    dist_max_pc: float = 100.0,
    tmag_max: float = 14.0,
    rad_max: float = 0.7,
    cache: Path | None = None,
) -> pd.DataFrame:
    """All-sky TIC query for nearby M dwarfs. ~90 s and ~90k rows with the defaults.

    astroquery.mast.Catalogs.query_criteria filter syntax:
      - a scalar selects equality:      lumclass="DWARF"
      - a 2-list selects a range:        Teff=[2300, 4000]
    Column names are the TIC names (case-sensitive): Teff, plx, Tmag, rad, mass, d ...
    Distance is not filterable directly with a range in every MAST version, so we
    filter on parallax (mas): d < 100 pc  <=>  plx > 10 mas.
    """
    if cache is not None and cache.exists():
        log.info("using cached TIC query %s", cache)
        return pd.read_parquet(cache)

    log.info("querying TIC (all-sky) ... this takes ~1-2 min")
    tab = Catalogs.query_criteria(  # type: ignore
        catalog="TIC",
        objType="STAR",
        lumclass="DWARF",
        Teff=[2300.0, teff_max],  # M dwarfs; TIC Teff is missing for many faint stars
        plx=[1000.0 / dist_max_pc, 1e4],  # mas
        Tmag=[0.0, tmag_max],
        rad=[0.0, rad_max],
    )
    df = tab[TIC_COLUMNS].to_pandas()
    df = df.rename(
        columns={
            "ID": "tic_id",
            "Tmag": "tmag",
            "Teff": "teff",
            "GAIA": "gaia_dr2_id",
            "d": "dist_pc",
        }
    )
    df["tic_id"] = df["tic_id"].astype("int64")
    log.info("TIC returned %d stars", len(df))
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache, index=False)
    return df


def stars_in_sector(df: pd.DataFrame, sector: int) -> pd.DataFrame:
    """Keep the rows that fall on a TESS camera/CCD in `sector` (geometric check via tess-point).

    tess_stars2px_function_entry returns one output row per (star, sector) combination
    over ALL sectors the star ever falls on (past and planned), so we filter on outSec.
    `trySector=sector` restricts the search to one sector and is much faster.
    """
    ids = df["tic_id"].to_numpy()
    out_id, _elon, _elat, out_sec, out_cam, out_ccd, out_col, out_row, _scinfo = (
        tess_stars2px_function_entry(
            ids, df["ra"].to_numpy(), df["dec"].to_numpy(), trySector=sector
        )
    )
    hit = pd.DataFrame(
        {
            "tic_id": out_id,
            "sector": out_sec,
            "camera": out_cam,
            "ccd": out_ccd,
            "col_pix": out_col,
            "row_pix": out_row,
        }
    )
    hit = hit[hit["sector"] == sector].drop_duplicates("tic_id")
    log.info("%d / %d stars fall in sector %d", len(hit), len(df), sector)
    return df.merge(hit, on="tic_id", how="inner").reset_index(drop=True)


def build_targets(sector: int, force: bool = False, **tic_kwargs) -> pd.DataFrame:
    out = DATA_DIR / f"sector_{sector}.parquet"
    if out.exists() and not force:
        log.info("%s exists, use --force to rebuild", out)
        return pd.read_parquet(out)
    allsky = query_tic_mdwarfs(
        cache=DATA_DIR / "tic_mdwarfs_allsky.parquet", **tic_kwargs
    )
    targets = stars_in_sector(allsky, sector)
    targets = targets.sort_values("tmag").reset_index(drop=True)
    out.parent.mkdir(parents=True, exist_ok=True)
    targets.to_parquet(out, index=False)
    log.info("wrote %s (%d stars)", out, len(targets))
    return targets


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    SECTOR = int(os.environ.get("SECTOR", 0))
    p.add_argument(
        "--sector",
        type=int,
        default=SECTOR,
        help="TESS sector number (retrieved from $SECTOR env var if not given)",
    )
    p.add_argument("--teff-max", type=float, default=4000.0, help="K")
    p.add_argument("--dist-max", type=float, default=100.0, help="pc")
    p.add_argument("--tmag-max", type=float, default=14.0, help="mag")
    p.add_argument("--rad-max", type=float, default=0.7, help="R_sun")
    p.add_argument("--force", action="store_true")
    a = p.parse_args()
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )
    t = build_targets(
        a.sector,
        force=a.force,
        teff_max=a.teff_max,
        dist_max_pc=a.dist_max,
        tmag_max=a.tmag_max,
        rad_max=a.rad_max,
    )
    print(
        t[
            [
                "tic_id",
                "ra",
                "dec",
                "tmag",
                "teff",
                "rad",
                "dist_pc",
                "contratio",
                "camera",
                "ccd",
            ]
        ].head(10)
    )
    print(
        f"\n{len(t)} targets; Tmag median {t.tmag.median():.1f}; Teff median {t.teff.median():.0f} K"
    )


if __name__ == "__main__":
    main()
