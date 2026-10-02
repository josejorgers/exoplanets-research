# Decisions and findings log

Short records of things we checked, what we found, and what it changed in the pipeline.
Newest first. Each entry says how to re-check it.

---

## 2026-10-02 — The 0.5-day biweight keeps ~90 % of short transits, ~75 % of 2.5-h ones

**Check** (`tests/test_detrend_inject.py`; plan §5 sanity check, not done before). A 2-mmag
box transit every 2 d, on a 1 %, 3-d sine trend with 1,000 ppm white noise, flattened with
`flatten_lc`. Fraction of the depth retained:

| duration | 0.5 h | 1.0 h | 1.5 h | 2.5 h |
|---|---|---|---|---|
| flat star, noise | 1.01 | 0.91 | 0.92 | 0.83 |
| 1 %/3 d trend, noise | 1.02 | 0.90 | 0.96 | 0.72 |

The biweight is robust but not blind: a transit filling ~20 % of the 12-h window pulls the
trend down a little. This matches the "window ≈ 3× the longest duration" rule (2.5 h × 3 =
7.5 h is close to the limit). A 1 %, **1-day** modulation is not removed at all (the residual
is anti-correlated with the box). That is by design, and the vetting variability test exists
for it.

**Consequences.** Depths of long transits (P ≳ 8 d around early M dwarfs) are underestimated
by up to ~25 %, and so are their `rp_earth` (by ~13 %). The injection-recovery run measures
what this costs in detections. A transit *fit* (follow-up) should be done on light curves
detrended with the transit masked.

## 2026-10-02 — QLP alternate-aperture fluxes cannot be used as a blend test

**Idea.** QLP files contain fluxes from three apertures: `DET_FLUX_SML`, `SAP_FLUX` (the
one our light curves use) and `DET_FLUX_LAG`. For a deblended light curve, a transit on the
target should have the same depth in every aperture. A depth that grows with the aperture
would then point to a neighbour.

**Check.** On deep eclipsing binaries that are clearly on their own, nearly uncontaminated
star (`contratio` < 0.2, SDE > 15, depth > 2 %), the depth ratio between the large and small
apertures ranged from **0.8 to 60**. One example: TIC 469939293 (contratio 0.001) gives
3,200 / 33,800 / 157,800 ppm. `DET_FLUX` agrees with `SAP_FLUX` (34,000 ppm), so QLP's own
detrending is not the cause. The alternate apertures are normalised differently, in a way the
files do not document.

**Decision.** No aperture test in `vet.py`. `notebooks/pipeline-6-known-walkthrough.ipynb`
§6 was corrected (its first version presented the WASP-145 aperture trend as a blend
test).

**Re-check.** Compare the depths at the TLS ephemeris in the three columns for QLP stars with
deep, isolated EBs. If a future QLP release documents a per-aperture dilution correction,
retry it.

## 2026-10-01 — Both light-curve products are deblended

**Finding.** TIC 381856446 lies 5″ from WASP-145 A (TOI-2206.01, a known hot Jupiter with a
10,110 ppm transit) and is 1.6 mag fainter. Its QLP `SAP_FLUX` shows the transit **34,300 ppm**
deep. A raw blend could show at most ~8,200 ppm. A deblended light curve, which subtracts the
TIC-predicted neighbour flux and renormalises, allows up to ~43,000 ppm. So QLP is
deblended, and SPOC PDCSAP is too (`CROWDSAP`).

**Consequences.**
- `known.py` uses the deblended form of the dilution argument: a neighbour needs
  $D_n \ge D_{obs} F_t / F_n$. The raw form, $D_{obs}(1 + F_t/F_n)$, reversed the
  conclusion for TOI-3282.01.
- `vet.py` no longer stores `rp_earth_undiluted` (it corrected twice). `rp_earth` is the
  physical radius.

## 2026-10-01 — Catalog ephemerides must be propagated with the catalog period

**Bug.** The first version of the ephemeris check compared epochs modulo *our* period. Over
the ~400–2,000 orbits since a catalogued epoch, a period difference of 0.002 d (one sector's
precision) accumulates to hours. Three confirmed planets (TOI-1080 b, TOI-244 b, TOI-4552 b)
were wrongly classified as `new_in_known_system`.

**Fix.** Predict with the catalog period, which is fitted over years: compare with events
spaced $\min(r, 1)\,P_{cat}$ apart, with $\sigma_t = \sqrt{\sigma_T^2 + (n\,\sigma_P)^2}$.

**How it was found.** Measuring the transit depth in our own light curve at both predicted
times showed that both ephemerides land on the same transits.

## 2026-10-02 — Scattered light: already removed by the quality flags

**Check.** In 20 random SPOC files, 54,564 cadences carry bit 13 (4096, scattered light) or
bit 15 (16384, bad calibration). PDCSAP is NaN on all of them, and none appear in
`data/lc_flat/`. For QLP, `detrend.py` keeps only `quality == 0`, which removes QLP's
scattered-light bit (2048).

**Decision.** No per-sector Data-Release-Notes masking in `config.toml`. The plan's §5 TODO
is closed.
