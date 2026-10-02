# Exoplanet transit search in TESS M dwarfs

A pipeline that searches one TESS sector for transiting planets around nearby M dwarfs
(< 100 pc), vets the signals, removes everything already known, and produces a ranked report
of new candidates. Milestone 1 of `docs/PLAN_milestone1.md`; first run on **sector 104** (status and results at
the top of that file: all 9 confirmed transiting planets in the target list found, 8 of them
passing vetting and labelled known; 19 new signals, none surviving follow-up).

## Pipeline

| step | module | what it does | output |
|---|---|---|---|
| 1 | `src/targets.py` | M dwarfs from the TIC that fall on TESS's detectors in the sector | `data/targets/sector_N.parquet` |
| 2 | `src/fetch.py` | download light curves (SPOC 2-min, else QLP FFI) | `data/lc_raw/`, `sector_N_manifest.parquet` |
| 3 | `src/detrend.py` | normalise, clip flares, biweight-flatten, trim orbit edges | `data/lc_flat/{tic}.parquet` |
| 4 | `src/search.py` | Transit Least Squares on every light curve, no threshold | `data/search/sector_N.parquet` |
| 5 | `src/vet.py` | 12 hard physical tests + 4 soft flags per signal | `data/search/sector_N_vetted.parquet` |
| 6 | `src/known.py` | cross-match with TOI, CTOI, confirmed planets, TESS & Gaia EBs, neighbours | `data/search/sector_N_known.parquet` |
| 7 | `src/report.py` | one PDF page per candidate + run summary | `data/reports/sector_N.pdf` |
| — | `src/inject.py` | injection-recovery: completeness of steps 3–5 | `data/injection/sector_N_injrec.parquet` |
| — | `src/followup.py` | other sectors, MCMC fit, difference image, Gaia, TRICERATOPS FPP | `data/followup/{tic}/` |

## Running

```bash
uv sync                                              # environment (Python 3.14)
echo "SECTOR=104" >> .env                            # default sector for every step
uv run python -m src.cli run --sector 104 --n-jobs 8 # the whole chain (hours: fetch and search dominate)
uv run python -m src.cli run --from vet              # re-run from a step
uv run python -m src.cli vet --rethreshold --snr-min 6   # any step with its own options
uv run python -m src.cli inject --n 1000 --n-jobs 10     # calibration (not part of `run`)
uv run --extra followup python -m src.followup --tics 306101698   # §12 follow-up of single candidates
uv run pytest -m "not slow"                          # unit tests (seconds)
uv run pytest -m slow                                # end-to-end on LHS 3844 (needs data/, ~1 min)
```

Every step skips work already on disk unless given `--force`. Logs go to
`data/logs/sector_N.log`, and package versions to `data/logs/sector_N_versions.json`.
Catalogs used by step 6 are kept as dated copies in `data/known/YYYY-MM-DD/`
(`--offline` reuses them).

Measured on sector 104 (12 cores): targets ~3 min, fetch ~1 h, detrend ~5 min,
search ~2–4 h, vet ~2.5 min, known ~10 s, report ~5 min (cold TLS cache),
1,000 injections ~2 h.

## Learning material

Each step has a walkthrough notebook that runs the method by hand on a few stars and explains
the parameters, the physics and the reasoning, for readers who know programming and math but
not astronomy:

| notebook | topic |
|---|---|
| `notebooks/pipeline-1-targets-walkthrough.ipynb` | the TIC, M dwarfs, which stars a sector observes |
| `notebooks/pipeline-2-fetch-walkthrough.ipynb` | TESS light-curve products and quality flags |
| `notebooks/pipeline-3-detrend-walkthrough.ipynb` | flattening without destroying transits |
| `notebooks/pipeline-4-search-walkthrough.ipynb` | Transit Least Squares, SDE, the periodogram |
| `notebooks/pipeline-5-vet-walkthrough.ipynb` | false positives and the physics of each vetting test |
| `notebooks/pipeline-6-known-walkthrough.ipynb` | catalogs, ephemeris propagation, blends and dilution |
| `notebooks/pipeline-7-report-walkthrough.ipynb` | how to read a candidate page; derived planet properties |
| `notebooks/pipeline-8-injection-recovery.ipynb` | completeness, and choosing the thresholds |
| `notebooks/pipeline-9-candidate-followup.ipynb` | from candidate to CTOI: other sectors, pixels, transit fit, false-positive probability |

Other documents:

- `docs/PLAN_milestone1.md`: the plan, with the reasoning behind each step.
- `docs/NOTES_decisions.md`: a log of what was checked and what it changed (with negative
  results).
- `docs/REVIEW_sector104.md`: the by-eye review of the sector-104 candidates.

## Layout

```
src/        pipeline modules (each runnable: python -m src.<step> --help)
tests/      pytest: synthetic-light-curve unit tests + an end-to-end test on LHS 3844
notebooks/  walkthroughs (see above)
docs/       plan, decision log, candidate review
data/       everything produced (git-ignored)
```
