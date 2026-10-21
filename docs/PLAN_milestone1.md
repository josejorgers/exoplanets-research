# Milestone 1 — Automated transit search for new candidates in one TESS sector

## Status (2026-10-19): milestone 1 run on sector 104

| plan step | module | status |
|---|---|---|
| 1 targets | `src/targets.py` | done: 4,710 M dwarfs |
| 2 fetch | `src/fetch.py` | done: 4,695 light curves (4,087 QLP, 608 SPOC) |
| 3 detrend | `src/detrend.py` | done; orbit-edge trim; scattered light handled by quality flags |
| 4 search | `src/search.py` | done: 4,693 TLS searches; second pass on the candidates |
| 5 vetting | `src/vet.py` | done: 13 hard tests + 4 soft flags (v3, after two reviews and injections) |
| 6 known objects | `src/known.py` | done: TOI, CTOI, planets, TESS & Gaia EBs, neighbours, siblings |
| 7 report | `src/report.py` | done: PDF + PNG per candidate + summary |
| 8 CLI & tests | `src/cli.py`, `tests/` | done: `python -m src.cli run`, 33 tests incl. end-to-end on LHS 3844 |
| §6 injection-recovery | `src/inject.py` | done: 1,000 injections; thresholds tuned (SDE ≥ 10, SNR ≥ 6.5) |
| §12 follow-up | `src/followup.py` | **automated** as a pipeline step (`followup`, notebook 10) on every new survivor: repetition (incl. 2P), null-calibrated difference image, fit, TRICERATOPS with background scenarios and Gaia-filled neighbours; verdict in the report |

**Definition of done** (§11):

- *full-sector run completes unattended*: yes (`python -m src.cli run --sector 104`);
- *recovers the known TOIs, flagged as known*: the target list contains 9 confirmed
  **transiting** planets. The search finds **all 9** at the right period, 8 pass vetting and are
  labelled `known_planet`, and TOI-122 b fails only the SNR threshold (SNR 5.6). TOI-6858.01 (PC)
  is found but fails SDE and SNR. (The other confirmed planets around these stars were found by
  radial velocities and do not transit; they are why some signals are classed
  `new_in_known_system`.)
- *ranked report of unknown signals*: yes, 13 new signals after vetting v3
  (`data/reports/sector_104.pdf`), graded in `docs/REVIEW_sector104.md`;
- *injection/recovery numbers documented*: yes (notebook 8, `docs/NOTES_decisions.md`). 50 %
  completeness at an expected SNR of 8.7 after vetting. The §6 target "> 80 % at 3 mmag, P < 10 d"
  is **not** met as written (73 %; ~65 % counting grazing orbits): half the losses are faint
  stars where 3 mmag is below SNR 7. At expected SNR ≥ 7 the same box reaches 85 %
  (notebook 8 §5.1–5.2).

**Science outcome:** one candidate, not yet submitted. Of the 13 "new" signals, the automatic
follow-up rejects 12 (7 eclipsing binaries or blends on neighbouring stars, 5 that do not
repeat, two of them systematics by eye). TIC 214311908 (P = 5.69 d, 4.9 R⊕) repeats over six
years, its bright neighbours are excluded by a PRF pixel fit, and its FPP is 0.09 ± 0.05
(0.45 ± 0.14 without the pixel exclusions): it passes the CTOI bar but not validation. An EB on
the 3.1″ neighbour (≈ 3.5σ) and a companion closer than ~1″ remain open
(`notebooks/candidate-TIC214311908.ipynb`, ephemeris for ground-based observers in
`data/followup/214311908/ephemeris.csv`). Next steps: decide on the CTOI, run further sectors.

Decision log: `docs/NOTES_decisions.md`. Walkthrough notebooks for every step: `notebooks/`.


## 0. Why this milestone

Candidate milestones were scored on *promise* (chance the output contains
something genuinely new) vs *difficulty* (data volume, tooling maturity,
domain knowledge needed):

| Milestone | Promise | Difficulty | Notes |
|---|---|---|---|
| **A. Nearby M dwarfs, one recent sector, ready-made FFI light curves, TLS search, cross-match against TOI/CTOI** | Medium–high | Low | Light curves already exist on MAST; M dwarfs give deep transits for small planets; recent sectors are less mined. **Chosen.** |
| B. Single-transit search in multi-sector stars | High | High | Needs custom detection logic, heavy vetting, long follow-up. |
| C. Faint stars (Tmag > 13.5) with own FFI photometry | High | High | Must do aperture photometry from TESSCut cutouts; systematics are hard. |
| D. Extra planets in known TOI systems | Medium | Medium | Small population; needs masking of known signals and careful residual analysis. |
| E. Kepler/K2 re-search | Low | Low | Thoroughly mined; unlikely to yield anything new. |

Milestone A is chosen. B/C/D are natural follow-ups that reuse the same pipeline skeleton.

**Deliverable:** a command-line pipeline that, for a given TESS sector, produces a
ranked PDF/HTML report of transit-like signals on nearby M dwarfs that are *not*
already listed as TOIs, CTOIs, confirmed planets, or known eclipsing binaries.

---

## 1. Concepts you need (read before coding)

Spend ~1 day here. Everything below is used in the pipeline.

- **Transit method**: a planet crossing its star dims it by (Rp/R*)². Earth-Sun ≈ 84 ppm; Earth-around-M-dwarf ≈ 1,000–3,000 ppm. That is why M dwarfs are the target.
- **Light curve**: flux vs time. TESS gives ~27 days per sector, one point every 2 min / 10 min / 200 s depending on product.
- **Detrending / flattening**: removing slow stellar variability and instrument drifts while preserving transit-length dips.
- **BLS / TLS**: algorithms that scan periods and look for periodic box- or transit-shaped dips. Output includes a period, depth, duration and a detection statistic (SDE — signal detection efficiency). SDE ≳ 7–9 is the usual "look at this" threshold.
- **False positives**: eclipsing binaries (EBs), blended EBs on nearby stars, stellar systematics, momentum dumps (TESS reaction-wheel desaturations cause spikes). Most detections are false positives; vetting is the bulk of the work.
- **Identifiers**: every TESS star has a **TIC ID**. This is the key for all cross-matches.

Reading (short list, ordered):
1. Lightkurve tutorial *"Identifying transiting exoplanet signals in a light curve"* — https://lightkurve.github.io/lightkurve/tutorials/3-science-examples/exoplanets-identifying-transiting-planet-signals.html
2. Lightkurve tutorial *"Removing noise with wotan-like flattening / Using the flatten() method"* — https://lightkurve.github.io/lightkurve/tutorials/
3. Winn (2010) *Transits and Occultations* — arXiv:1001.2010 (sections 1–3 are enough)
4. Hippke & Heller (2019) *Transit Least Squares* — arXiv:1901.02015
5. Hippke et al. (2019) *wotan* — arXiv:1906.00966
6. TESS Data Release Notes for the sector you choose — https://archive.stsci.edu/tess/tess_drn.html

---

## 2. Environment setup

Project already uses `uv`, Python ≥ 3.13, lightkurve, astropy, numpy, pandas, matplotlib.

Add dependencies:

| Package | Purpose |
|---|---|
| `astroquery` | Query MAST (TIC), Vizier, NASA Exoplanet Archive |
| `wotan` | Detrending (biweight filter) |
| `transitleastsquares` | TLS transit search |
| `tqdm` | Progress bars |
| `pyarrow` | Fast Parquet I/O for tables |
| `joblib` | Simple parallelism over stars |

Optional (later milestones): `triceratops`, `batman-package`, `exoplanet`/`juliet`, `eleanor`.

Check Python 3.13 compatibility of `transitleastsquares` and `wotan`; if a wheel fails to build, pin Python to 3.12 in `.python-version` — both packages are pure Python/numba and normally fine.

Create a `.env` or `config.toml` with: data directory, sector number, magnitude cut, TLS thresholds, number of workers.

Proposed layout:

```
exoplanets-research/
  config.toml
  src/
    __init__.py
    targets.py      # step 3
    fetch.py        # step 4
    detrend.py      # step 5
    search.py       # step 6
    vet.py          # step 7
    known.py        # step 8
    report.py       # step 9
    cli.py          # orchestrates: python -m src.cli run --sector 104
  data/
    targets/        # target lists (parquet)
    lc_raw/         # downloaded light curves (fits)
    lc_flat/        # detrended light curves (parquet)
    search/         # TLS results (parquet, one row per star)
    known/          # TOI, CTOI, EB catalogs
    reports/
  tests/
  notebooks/        # exploratory only
```

Each step reads from and writes to `data/` so any step can be re-run in isolation.

---

## 3. Step 1 — Build the target list (`src/targets.py`)

Goal: `data/targets/sector_N.parquet`, one row per nearby M dwarf that TESS
observed in sector N, with the stellar parameters later steps need.

Reference implementation: `src/targets.py` (run with
`uv run python -m src.targets --sector 104`). Verified 2026-09-13: 2 min 42 s
end-to-end, 4,710 stars for sector 104. Everything below explains that file
line by line so you can modify it.

### 3.1 Choose the sector

Pick the most recently released sector that has QLP or TESS-SPOC light curves
on MAST (see below). Record the sector dates and known issues
from its Data Release Notes. The chosen sector goes in `.env` as `SECTOR=104`.

#### How to find the most recent usable sector

Different products lag the observations by different amounts, so "latest sector" depends on which light-curve product you use. Check, in order:

1. **What has been observed** — TESS observation schedule: https://tess.mit.edu/observations/ (current sector and dates).
2. **What has Data Release Notes** (i.e. is released on MAST): https://archive.stsci.edu/tess/tess_drn.html — the largest sector number listed.
3. **Which product exists for that sector** — each HLSP stores files in per-sector folders named `s00NN`. Open the folder listing and look at the largest `s00NN`:
   - QLP: https://archive.stsci.edu/hlsps/qlp/
   - TESS-SPOC: https://archive.stsci.edu/hlsps/tess-spoc/
   - TGLC: https://archive.stsci.edu/hlsps/tglc/
   - SPOC 2-min (not HLSP): https://archive.stsci.edu/tess/bulk_downloads/bulk_downloads_ffi-tp-lc-dv.html — largest `tesscurl_sector_NN_lc.sh`.
   From a terminal: `curl -sL https://archive.stsci.edu/hlsps/qlp/ | grep -oE 'href="s[0-9]{4}/"' | sort -u | tail -3` — use four digits (`s0104`), not `s00NN`, or sectors ≥ 100 are missed. Ignore folders that are not plausible sector numbers (e.g. `s1751` in tess-spoc is a multi-sector/other product, not a sector).
4. **Confirm with lightkurve** on a star in the continuous viewing zone (observed in most sectors), e.g. TOI-700 (`TIC 150428135`, south) or TIC 261136679 (north): `search_lightcurve("TIC 150428135").table` → group by `author`, take max `sequence_number` (= sector). Note the star is only in one hemisphere's CVZ, so it lags by up to a year when TESS points the other way — step 3 is the authoritative check.
5. Pick the sector, then read its Data Release Notes before using it.

Steps 2–4 are automated in `scripts/latest_sector.sh` (`-p qlp|tess-spoc|tglc`). It also estimates whether the newest folder is still being uploaded by counting TIC-prefix folders and comparing with the previous sectors (ratio ≥ 0.85 → complete). Caveat: the number of stars per sector varies a lot with the field (galactic plane sectors have ~2× more), so treat a ratio of 0.5–0.85 as "check by hand", and always require the DRN to be published. Run it once, wait ~5 minutes (it walks thousands of folder listings).

Snapshot on 2026-09-13: DRN up to sector 104; SPOC 2-min light curves up to sector 106; QLP folder exists for sector 105 but is **still uploading** (3,195 TIC-prefix folders vs 11,748 for s104 and 21,874 for s099) and has no DRN yet; TESS-SPOC up to sector 85; TGLC up to sector 55. Recommended for milestone 1: **sector 104 with QLP** (re-run the script later and switch to 105 once complete; use 85 if TESS-SPOC quality is preferred). Typical lags: SPOC ≈ 1 month, QLP ≈ 1–2 months, TESS-SPOC ≈ 1–1.5 years.

### 3.2 Query the TIC with astroquery — the API, precisely

**Which class.** The TIC is a *catalog*, so it is queried through
`astroquery.mast.Catalogs`. `astroquery.mast.Observations` is for
*observation products* (FITS files) and does not accept `catalog="TIC"`.

**Import gotcha.** `astroquery` does *not* import its sub-packages on
`import astroquery`; `astroquery.mast` is only defined after you import it
explicitly. This is why the notebook cell fails with
`AttributeError: module 'astroquery' has no attribute 'mast'`. Always write:

```python
from astroquery.mast import Catalogs      # correct
# import astroquery; astroquery.mast...  # AttributeError
```

**Filter syntax of `Catalogs.query_criteria`.** Every keyword argument that is
not `catalog`/`coordinates`/`radius` is treated as a column filter on the TIC:

| You write | Meaning |
|---|---|
| `lumclass="DWARF"` | equality |
| `Teff=[2300, 4000]` | inclusive range `2300 ≤ Teff ≤ 4000` (a 2-element list) |
| `objType="STAR"` | equality (excludes `EXTENDED` galaxies) |

Column names are the TIC's own, case-sensitive (`Teff`, `Tmag`, `plx`, `rad`,
`mass`, `d`, `contratio`, `lumclass`, `objType`, `GAIA`, ...). The full list is
`Catalogs.query_criteria(catalog="TIC", ID=150428135).colnames` or Table 1 of
Stassun et al. 2019 (arXiv:1905.10694). A wrong column name returns an HTTP
error from MAST, not a Python error.

**Recommended query (all-sky):**

```python
from astroquery.mast import Catalogs

tab = Catalogs.query_criteria(
    catalog="TIC",
    objType="STAR",
    lumclass="DWARF",
    Teff=[2300.0, 4000.0],   # M dwarfs (K7 boundary ~4000 K)
    plx=[10.0, 1e4],         # parallax in mas: plx > 10 mas  <=>  d < 100 pc
    Tmag=[0.0, 14.0],        # QLP/TESS-SPOC go to ~13.5; margin for later cuts
    rad=[0.0, 0.7],          # R_sun; removes giants misclassified as DWARF
)
```

Notes on the cuts:
- Distance is filtered through **parallax** (`plx`), not `d`; `plx` is always
  populated where `d` is, and the conversion is `d[pc] = 1000 / plx[mas]`.
- No coordinate filter is applied. Restricting to the sector's sky footprint
  with `coordinates=`/`radius=` is *not* practical (a sector is a 24°×96° strip,
  not a cone), so query all-sky once and cache it — the sector filter is done
  locally in §3.3.
- The result is an `astropy.table.Table`. Measured: **88,934 rows, ~90 s**.
  Convert with `tab[cols].to_pandas()` and cache to
  `data/targets/tic_mdwarfs_allsky.parquet` so the query runs once per
  parameter set, not once per sector.
- `ID` (the TIC ID) comes back as a string in some astroquery versions —
  cast to `int64` (`df["tic_id"] = df["ID"].astype("int64")`).

Columns to keep (renamed in the parquet): `ID→tic_id, ra, dec, Tmag→tmag,
Teff→teff, rad, mass, d→dist_pc, plx, contratio, lumclass, objType,
GAIA→gaia_dr2_id, pmRA, pmDEC`.

**Failure modes.** MAST times out occasionally (`TimeoutError` /
`RemoteServiceError`) — rerun; the query is idempotent. If it happens
repeatedly, split the query into two `Tmag` ranges (e.g. `[0,12]` and
`[12,14]`) and concatenate.

### 3.3 Filter to the stars TESS actually observed in the sector

The TIC has no notion of sectors. Use **tess-point** (package `tess-point`,
module `tess_stars2px`), which applies the spacecraft pointing model and
returns camera/CCD/pixel for every (star, sector) pair. Signature:

```python
from tess_stars2px import tess_stars2px_function_entry

(out_id, out_elon, out_elat, out_sec, out_cam, out_ccd,
 out_col, out_row, scinfo) = tess_stars2px_function_entry(
    tic_ids,            # np.ndarray of int64
    ra_deg, dec_deg,    # np.ndarray, degrees (ICRS)
    trySector=104,      # restrict to one sector: much faster than all sectors
)
```

- Inputs are arrays: one call for all ~89k stars (**~70 s** for sector 104).
- Output arrays have one row per *(star, sector)* hit, and stars that do not
  fall on any CCD are simply absent. Filter `out_sec == sector`, then
  `drop_duplicates("tic_id")` (a star at a CCD edge can be reported twice).
- tess-point's built-in pointing table currently goes to sector 134
  (`tess_stars2px.max_sector`); if you ever pick a sector beyond that,
  `uv lock --upgrade-package tess-point`.
- Do not confuse `tess_stars2px_function_entry` (Python API) with
  `tess_stars2px.py` used as a CLI (`python -m tess_stars2px -c ra dec`) — the
  CLI is fine for one star, useless for thousands.

Alternative (slower, but ground truth): `lightkurve.search_lightcurve(f"TIC {tic}", sector=N)`
per star and keep those with a result. ~1–2 s per star, so use it only for
spot checks (§3.5) or as a fallback in step 2 (stars with no light curve are
dropped there anyway).

### 3.4 Save

`data/targets/sector_N.parquet`, sorted by `tmag`, columns:
`tic_id, ra, dec, tmag, teff, rad, mass, dist_pc, plx, contratio, lumclass,
objType, gaia_dr2_id, pmRA, pmDEC, sector, camera, ccd, col_pix, row_pix`.

Also keep `data/targets/tic_mdwarfs_allsky.parquet` (the cached §3.2 result;
delete it to force a fresh TIC query when you change the cuts).

Measured result, sector 104: **4,710 stars** (Tmag median 12.8, Teff median
3,450 K; ~1,150 per camera). Tighter cuts: `Tmag < 13.5 → 3,675`,
`Tmag < 13 → 2,647`. Start the pipeline with `--tmag-max 13`; widen later.

### 3.5 Sanity checks (do all of them once)

1. `df.teff.isna().sum() == 0` and `df.camera.value_counts()` roughly equal
   across cameras 1–4 (a strongly unequal split means a wrong sector or a
   pointing-table problem).
2. Histograms of `tmag` and `teff`; `teff` should pile up at 3,200–3,600 K.
3. For 3 random `tic_id`, `lightkurve.search_lightcurve(f"TIC {tic}", sector=N)`
   must return at least one row (author `QLP`, `exptime` 200 s for sector 104).
   Verified for TIC 201255893, 360876935, 179823137.
4. Open one on ExoFOP (`https://exofop.ipac.caltech.edu/tess/target.php?id=<tic>`)
   and confirm sector N appears in its "TESS Sectors" line.
5. Regression: a known TOI in that sector must be in the list (find one with
   the ExoFOP TOI table filtered on sector; e.g. for a southern sector, check
   whether TOI-700 / TIC 150428135 was observed via tess-point).

Sources: TIC v8.2 description — Stassun et al. 2019, arXiv:1905.10694. Gaia DR3 (for later, better stellar parameters) — https://gea.esac.esa.int/archive/.

---

## 4. Step 2 — Download light curves (`fetch.py`)

Goal: one light-curve file per star, cached locally.

1. For each TIC ID, call `lightkurve.search_lightcurve(f"TIC {tic}", sector=N, author=["SPOC", "TESS-SPOC", "QLP"])`. Preference order: **SPOC** (2-min, best quality) → **TESS-SPOC** (FFI, PDCSAP flux, good systematics correction) → **QLP** (FFI, KSPSAP flux). Pick the first available.
2. Download with `.download()` to `data/lc_raw/`. Lightkurve caches in `~/.lightkurve/cache`; set `LIGHTKURVE_CACHE_DIR` env var or symlink to keep it inside the project.
3. Use the `quality` bitmask (`quality_bitmask="default"`) to drop flagged cadences. Also remove NaNs.
4. Record in a manifest (`data/targets/sector_N_manifest.parquet`): `tic_id, author, cadence_s, n_points, file_path, status` (ok / not_found / error).
5. Run with a small `joblib` pool (≤ 4 workers; MAST rate-limits) and retry on transient HTTP errors.

Alternative for bulk: MAST publishes per-sector bulk download shell scripts for SPOC and TESS-SPOC (https://archive.stsci.edu/tess/bulk_downloads.html). Faster for whole sectors, then filter by TIC ID locally.

Sanity checks: plot 5 random raw light curves; the manifest should show ≥ 80 % `ok`.

---

## 5. Step 3 — Detrend (`detrend.py`)

Goal: a flat light curve with transits preserved.

1. Normalize flux to 1 (`lc.normalize()`).
2. Sigma-clip *upward* outliers only (flares, cosmic rays): e.g. `remove_outliers(sigma_upper=4, sigma_lower=inf)`. Never clip downward — you would clip transits.
3. Flatten with **wotan** `flatten(time, flux, method="biweight", window_length=0.5, return_trend=True)`. Rule of thumb: window ≈ 3× the longest transit duration you want to keep; 0.5 d is a good default for M dwarfs (durations ~1–2 h). Also try `lightkurve.flatten(window_length=...)` as a cross-check.
4. Mask known bad time ranges from the Data Release Notes (scattered light, momentum dumps) — store them in `config.toml` per sector.
5. Save `data/lc_flat/{tic}.parquet` with `time, flux, flux_err, trend`.

Sanity checks: standard deviation of flattened flux should be close to the median `flux_err`; inject a fake 2 mmag box dip in one star and confirm it survives flattening.

**Orbit-edge fix (found in `notebooks/pipeline-4-search-walkthrough.ipynb`, §4/§8; implemented 2026-09-20).**
The raw QLP/SPOC flux climbs by ~10 % over the first day of each orbit (scattered
light as the spacecraft comes off perigee). The biweight window is one-sided at
the edge and leaves a ~2 % residual dip there, which TLS then reports as a
high-SDE "transit" (e.g. TIC 6283622: SDE 23–28, P = `period_max`, 1 of 2
transits with data; same for LHS 3844 in the second pass). This was 2 of the 5
stars searched by hand, so it was the most common false positive in the sector.

What `detrend.py` now does (`edge_mask`, `--edge-trim`, default 0.25 d):
split the light curve at gaps > 0.5 d (the orbits), flatten on *all* points,
then drop the first/last `edge_trim` days of each orbit. The count goes to
`n_trimmed` in the detrend summary. Measured on TIC 6283622 (deepest 1.7 h
running-median dip near the edges, vs 0.30 % elsewhere): untrimmed 2.3 %;
trim 0.25 d *before* flattening 1.9 % (the ramp just moves the edge — do not
do this); trim 0.25 d *after* flattening 0.36 % ≈ noise. `search.py` on the
re-detrended star: SDE 27.8 → 6.6 (noise). `edge_trim ≈ window_length / 2`
is the region where the window is one-sided; scale it if you change the window.
Cost: ~360 of 8,800 points per star (4 %).

Scattered-light windows (item 4 above): **resolved without extra code** (checked
2026-10-02). SPOC sets PDCSAP to NaN on every cadence flagged with bit 13 (4096,
scattered light) or bit 15 (16384, bad calibration), and `detrend.py` keeps only
`quality == 0` for QLP. In 20 SPOC files, 0 of 54,564 flagged cadences survive into
`data/lc_flat/`. See `docs/NOTES_decisions.md`. Whatever survives still needs the
vetting rules in §7.

Source: wotan docs — https://github.com/hippke/wotan

---

## 6. Step 4 — Transit search (`search.py`)

Goal: one row per star with the best TLS signal.

1. For each flattened light curve run `transitleastsquares.transitleastsquares(time, flux, flux_err).power(...)` with:
   - `R_star`, `M_star` from the TIC (improves duration prior)
   - `period_min = 0.5 d`, `period_max = (time span)/2` (need ≥ 2 transits; ≈ 13 d for one sector)
   - `oversampling_factor=3`, `duration_grid_step=1.1`
   - `use_threads` = 1 per star, parallelize across stars with `joblib`
2. Record: `period, T0, depth, duration, SDE, snr, odd_even_mismatch, transit_count, distinct_transit_count, per-transit depths, FAP`.
3. Optional second pass: `--second-pass` reads the first-pass table, masks each star's best signal (`transit_mask`, 2× duration) and re-runs TLS to find a second planet. Output has the same columns, in `data/search/pass2/sector_N_pass2.parquet`.
   Use `period_min = 0.2 d`, not 0.5: ultra-short-period planets around M dwarfs are a real population (LHS 3844 b at 0.463 d is in sector 104) and short periods add few grid points. Do **not** fix the orbit-edge false positive (§5 TODO) in this step; only store the fields vetting needs to reject it (`distinct_transit_count`, `empty_transit_count`, `period_unc`, `period_max`, `transit_depths`).
4. Save all rows to `data/search/sector_N.parquet` (don't threshold here — thresholds are applied in vetting so you can tune them).

Runtime: TLS on a 27-day, 10-min-cadence light curve is ~10–60 s per star on one core. 1,000 stars × 30 s / 8 cores ≈ 1 h. Start with 50 stars to time it.

Sanity check — **injection/recovery** (mandatory before believing anything): inject synthetic transits (`batman` or a simple box) of known period/depth into ~50 real light curves, run steps 3–4, and measure the fraction recovered as a function of depth and period. Aim for > 80 % recovery at 3 mmag / P < 10 d. This tells you what your pipeline can and cannot see.

Source: TLS docs — https://transitleastsquares.readthedocs.io/

---

## 7. Step 5 — Automated vetting (`vet.py`)

Goal: reduce hundreds of "signals" to a handful worth human time. Implement each test as a function returning a pass/fail plus a value, and store all of them (no early exit).

| Test | Rule of thumb | Catches |
|---|---|---|
| SDE threshold | SDE ≥ 9 (tune with injection tests) | noise |
| SNR | snr ≥ 7 | noise |
| Transit count | ≥ 2 distinct transits with individually significant depth | single outliers |
| Grid-edge period | `period` not within a few grid steps of `period_max` (stored per star by `search.py`) and `period_unc` finite; with `distinct_transit_count == 1` this is the orbit-edge artefact of §5 — see the walkthrough notebook §4 | start-of-orbit scattered-light residual, single events |
| Odd/even depth | odd and even depths agree within 3σ | EBs at 2× period |
| Secondary eclipse | no significant dip at phase 0.5 | EBs |
| Duration vs stellar density | duration consistent with circular orbit around a star of TIC radius/mass (within ×3) | blends, wrong host |
| Depth | depth < 5 % (larger implies stellar companion for M dwarfs) | EBs |
| Transit shape | TLS `transit_shape`/V-shape flag; or compare box vs trapezoid fit | grazing EBs |
| Momentum-dump proximity | no transit within ±0.1 d of a momentum dump | systematics |
| Centroid (optional, needs TPF) | flux-weighted centroid does not shift during transit | blends |
| Contamination | TIC `contratio` < 0.5 | nearby bright neighbours |

Output: `data/search/sector_N_vetted.parquet` with all metrics and a boolean `passes_all`. Also compute a simple score (e.g. SDE × passes) for ranking.

Sources (the definitive vetting logic, adapt rather than reinvent):
- Thompson et al. 2018, Kepler Robovetter — arXiv:1710.06758
- Coughlin et al. 2016 — arXiv:1512.06149
- Twicken et al. 2018, Data Validation — arXiv:1803.04526
- LEO-vetter — https://github.com/mkunimoto/LEO-vetter (Python implementation of many of these tests; consider using it directly)
- DAVE — https://github.com/exoplanetvetting/DAVE
- TRICERATOPS (false-positive probability, later milestone) — Giacalone et al. 2021, arXiv:2002.00691, https://github.com/stevengiacalone/triceratops

---

## 8. Step 6 — Remove already-known objects (`known.py`)

Goal: keep only signals that are *not* already catalogued.

1. **TOI table** — download CSV from ExoFOP: https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv
2. **CTOI table** — https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv
3. **Confirmed planets & TOIs** — NASA Exoplanet Archive TAP: https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html (tables `ps`, `toi`); or `astroquery.ipac.nexsci.nasa_exoplanet_archive`.
4. **TESS Eclipsing Binary catalog** — Prša et al. 2022, arXiv:2110.01620; Vizier `J/ApJS/258/16`.
5. Cross-match by **TIC ID** first. Then, for candidates whose TIC ID is not listed, cross-match by *ephemeris*: same period (or ½, 2×) and a T0 that phases up, against known objects within ~1 arcmin (catches blends from a neighbouring known EB). Use `astropy.coordinates.SkyCoord.search_around_sky`.
6. Flag each candidate: `known_toi`, `known_ctoi`, `known_planet`, `known_eb`, `neighbour_match`. Keep everything in the table; the report just sorts known ones to the bottom.

Re-download these catalogs every run (they change weekly); keep dated copies in `data/known/`.

---

## 9. Step 7 — Report (`report.py`)

Goal: a document you can scan in 10 minutes.

For each candidate that passes vetting and is not known, one page with:
- Header: TIC ID, Tmag, Teff, R*, distance, sector, light-curve author.
- Panel 1: raw light curve with trend overlaid.
- Panel 2: flattened light curve with transit times marked.
- Panel 3: phase-folded light curve with binned points and TLS model.
- Panel 4: odd vs even folded transits.
- Panel 5: phase 0.5 (secondary) zoom.
- Panel 6: TLS periodogram (SDE vs period).
- Table: all vetting metrics, planet radius estimate (`sqrt(depth) × R*`), equilibrium temperature estimate, links to ExoFOP and the TIC page.

Output: `data/reports/sector_N.pdf` (matplotlib `PdfPages`) plus a summary CSV. Sort by score.

Also write a run summary: number of targets, downloaded, searched, passing each vetting stage, known, new.

---

## 10. Step 8 — CLI and reproducibility (`cli.py`)

- `exo targets --sector N`, `exo fetch`, `exo detrend`, `exo search`, `exo vet`, `exo known`, `exo report`, and `exo run --sector N` for the full chain.
- Each step skips stars whose output already exists unless `--force`.
- Log to `data/logs/sector_N.log`.
- Pin package versions in `uv.lock`; record the versions and config in the report header.

Tests (`tests/`):
- Unit tests for each vetting function using synthetic light curves.
- An end-to-end test on 3 stars including one with a known TOI (it must be detected *and* flagged as known) — good known test targets: TOI-700 (TIC 150428135), LHS 3844 (TIC 410153553), TOI-270 (TIC 259377017). Search these in a sector where they were observed.

---

## 11. Suggested schedule

| Week | Work |
|---|---|
| 1 | Reading (§1), environment (§2), target list (§3). Verify 3 known TOIs are in the target-list logic when their sector is chosen. |
| 2 | Fetch (§4) and detrend (§5). Plots of 20 stars reviewed by hand. |
| 3 | Search (§6) on 50 stars, then injection/recovery. Tune window length and SDE threshold from the results. |
| 4 | Vetting (§7) and known-object cross-match (§8). Confirm known TOIs get recovered and flagged. |
| 5 | Report (§9), CLI (§10), full-sector run, hand review. |

Definition of done: a full-sector run completes unattended, recovers the known TOIs in that sector (flagged as known), and produces a ranked report of unknown signals with injection/recovery numbers documented.

---

## 12. What happens after a candidate looks real

1. Check every sector the star has been observed in (`search_lightcurve` without `sector=`) and confirm the signal repeats.
2. Look at the target pixel file and centroids; check nearby Gaia sources within the aperture.
3. Fit a transit model (`batman` + `emcee`, or `juliet`) for radius, impact parameter, and density check.
4. Run TRICERATOPS for a false-positive probability.
5. If it holds up, publish it (a refereed paper, or an AAS Research Note citing published methods), then upload it as a **CTOI** on ExoFOP with the paper link, ephemeris and plots. ExoFOP only accepts published candidates and gives upload accounts to contributing astronomers (https://exofop.ipac.caltech.edu/tess/candidate_help.php).

---

## 13. Reference list (for digging deeper)

**Data & catalogs**
- MAST portal — https://mast.stsci.edu/
- TESS bulk downloads — https://archive.stsci.edu/tess/bulk_downloads.html
- TESS Science Support Center — https://heasarc.gsfc.nasa.gov/docs/tess/
- TESS Data Release Notes — https://archive.stsci.edu/tess/tess_drn.html
- TESSCut (FFI cutouts) — https://mast.stsci.edu/tesscut/
- HLSP light curves: QLP (Huang et al. 2020, arXiv:2011.06459), TESS-SPOC (Caldwell et al. 2020, arXiv:2010.00019), TGLC (Han & Brandt 2023, arXiv:2301.03704), eleanor (Feinstein et al. 2019, arXiv:1903.09152), DIAmante (Montalto et al. 2020, arXiv:2009.05946)
- TIC v8.2 — Stassun et al. 2019, arXiv:1905.10694
- Gaia DR3 archive — https://gea.esac.esa.int/archive/
- `tess-point` — https://github.com/christopherburke/tess-point

**Known objects**
- NASA Exoplanet Archive — https://exoplanetarchive.ipac.caltech.edu/
- ExoFOP-TESS (TOI/CTOI) — https://exofop.ipac.caltech.edu/tess/
- TESS EB catalog — Prša et al. 2022, arXiv:2110.01620
- Planet Hunters TESS — https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess; Eisner et al. 2021, arXiv:2011.13944

**Methods**
- BLS — Kovács, Zucker & Mazeh 2002, arXiv:astro-ph/0206099
- TLS — Hippke & Heller 2019, arXiv:1901.02015; https://transitleastsquares.readthedocs.io/
- wotan — Hippke et al. 2019, arXiv:1906.00966; https://github.com/hippke/wotan
- Winn 2010, *Transits and Occultations* — arXiv:1001.2010
- SPOC pipeline — Jenkins et al. 2016, SPIE 9913
- Faint-star QLP search — Kunimoto et al. 2022, arXiv:2202.03656
- Single-transit searches — Gill et al. 2020, arXiv:2005.00006; Villanueva et al. 2019, arXiv:1811.09757

**Vetting & validation**
- Kepler Robovetter — Thompson et al. 2018, arXiv:1710.06758; Coughlin et al. 2016, arXiv:1512.06149
- Data Validation — Twicken et al. 2018, arXiv:1803.04526
- LEO-vetter — https://github.com/mkunimoto/LEO-vetter
- DAVE — https://github.com/exoplanetvetting/DAVE
- TRICERATOPS — Giacalone et al. 2021, arXiv:2002.00691
- Astronet / exoplanet-ml — Shallue & Vanderburg 2018, arXiv:1712.05044; https://github.com/google-research/exoplanet-ml

**Modelling (post-candidate)**
- batman — https://lkreidberg.github.io/batman/
- exoplanet — https://docs.exoplanet.codes/
- juliet — https://juliet.readthedocs.io/

**Tutorials**
- Lightkurve tutorials — https://lightkurve.github.io/lightkurve/tutorials/index.html
- MAST notebooks (TESS section) — https://github.com/spacetelescope/notebooks
- Astroquery MAST docs — https://astroquery.readthedocs.io/en/latest/mast/mast.html
