# Decisions and findings log

Short records of things we checked, what we found, and what it changed in the pipeline.
Newest first. Each entry says how to re-check it.

---

## 2026-10-19 — Follow-up in the pipeline; honest FPPs; completeness explained

**1. Follow-up is a pipeline stage.** `python -m src.cli run` now runs `followup` between
`known` and `report` on every survivor with `known_class` in {new, new_in_known_system}:
repetition (all sectors, 2P search), difference image, and for the ones still standing the fit
and TRICERATOPS. `verdict()` in `src/followup.py` is shared by notebook 9, notebook 10 and the
report. Slow steps are reused when P and T0 are unchanged; `--verdict-only` recomputes verdicts
from stored summaries. On sector 104 it reproduces every by-eye rejection (12 of 12). Re-check:
`uv run --extra followup python -m src.followup --sector 104 --verdict-only`.

**2. Difference image: noise from pseudo-transits.** The old per-transit scatter gave S/N 7–10
on pure noise with 4 transits. The noise now comes from up to 60 difference images at the same
ephemeris shifted off-transit; a loss counts above max(5, the null's 99th percentile) anywhere
or max(3, …) within 3 px; otherwise `status = inconclusive`. Changed: 372868160 off → inconclusive
(verdict unchanged: it does not repeat); 214311908 S/N 8.0 → 4.3 (threshold 4.25). Tests:
`tests/test_diffimage.py`.

**3. TRICERATOPS: background scenarios and neighbour parameters.** TRILEGAL had failed because
its server omits the ZeroSSL intermediate certificate; `src/background.py` fetches it once,
keeps verification on, and caches each population. TRICERATOPS gives neighbours without TIC
parameters solar values; `src/neighbours.py` fills them from Gaia DR3 (the red giant near TIC
214311908 went from 1 to ~93 R☉), drops TIC duplicates, and applies optional per-target
exclusions from a PRF pixel fit (`data/followup/<tic>/pixel_excluded.json`). FPPs are the mean
± std of 5 runs (scatter up to ±0.14). A verdict never says "candidate" without background
scenarios. TIC 214311908: 0.64 → 0.45 ± 0.14 (no exclusions) → 0.09 ± 0.05 (pixel exclusions);
control LHS 3844 b: 0.015 (one lucky run) → 0.07 ± 0.05. Tests: `tests/test_neighbours.py`.

**4. Why the plan's completeness goal reads 73 % (notebook 8 §5.1–5.2).** The goal
"> 80 % at ≥ 3 mmag, P < 10 d" fixes depth, not noise. Of 362 such injections, 98 are lost: 53
below expected SNR 7 (faint stars, median scatter 23,000 ppm per point), 22 search losses (12
where the star's own stronger signal wins; the second pass rescues 1), 23 vetting losses (14 just
under the thresholds). The same box reaches 85 % at expected SNR ≥ 7 and 91 % at ≥ 10. Grazing
orbits (b > 0.9, 14 % of transiting orbits), left out of the main run, pass 7 % of the time;
excluding them overstates completeness by ~4 points overall and ~9 in the goal box. Better goal
for later milestones: ≥ 80 % at expected SNR ≥ 7, P < 10 d.

**5. Short gaps: no ramps, one QLP artefact.** Stacking 2,000 stars around the five 0.1–0.26-d
gaps shows ≤ 170 ppm residual ramps (negligible). Every QLP light curve has a ~1,100-ppm, 5-h
dip at BTJD 4197.62–4197.84 (in the raw QLP flux, not in SPOC); it attracts 1.4× more weak
TLS peaks but no signal with SDE ≥ 10. No mask (notebook 3 §8.3).

---

## 2026-10-04 — Vetting v3 and a broken follow-up repeat check

A review of the primer article against the code found two vetting cuts that cost planets or
let blends through, and a silent bug in the follow-up. Each change is tested in
`tests/test_vet.py` / `tests/test_followup.py`.

**1. Follow-up "narrow-window" searches were mostly full searches (bug).** TLS replaces its
period grid by the default one (R = M = 1, *no period limits*) when fewer than 100 grid periods
fall inside [period_min, period_max], and `followup._tls_window` suppressed the warning. A ±2 %
window around P = 4.5 d holds ~25 periods in one sector, so every per-sector search of TIC
394585043, 229804549 and 306101698, and two of 304448456, searched all periods (the stored
results show periods far outside the window, e.g. 1.13 d for a 4.50 d signal). Combined
multi-sector searches were not affected (long baseline → dense grid). **Fix:** raise
`oversampling_factor` until the window holds ≥ 100 periods (`window_oversampling`), and raise
an error if TLS still leaves the window.

**2. Repetition had no notion of a sector's sensitivity.** A 30-minute sector with 1,100
points cannot reach SNR 5 for a signal found at SNR 7–8, yet counted as "does not repeat".
**Fix:** `sector_search` stores `snr_expected` (discovery depth / sector scatter ×
√(in-transit points)); only sectors with `snr_expected` ≥ 7.5 can refute, a detection anywhere
counts, and `repetition()` returns `testable = False` when no sector is informative.

**2b. Eclipsing binaries at half their period need their own search.** A ±2 % window at P can
never return 2P, so the old "returns 2P" count only worked *because of* bug 1 (TIC 304448456's
other sectors were full searches that peaked at 2P). **Fix:** `alias_search` repeats the narrow
search around 2P per sector; a sector shows unequal eclipses when its 2P depth exceeds the P
depth by ≥ 20 % and ≥ 2σ (the ratio alone fired on half the sectors of TIC 306101698, which
repeats at P: TLS's best depth at 2P is biased high at modest SNR).

**2c. Off target needs the catalog position too.** The difference-image offset was measured
from the out-of-transit light centroid only, which neighbours pull away from the star. TIC
214311908 (contratio 2.2): 0.76 px from the light centroid, 0.30 px from the catalog position.
Now both must be > 0.5 px. No other follow-up verdict changed.

**3. Duration: no lower bound, upper bound 1.5 (was 1/3–3).** 8 of 412 correctly found
injections had TLS durations 0.09–0.32 of the expected one (TLS underestimates short
durations; most had small b) and were rejected. On the other side, 99 % of found injections
are ≤ 1.01, the 8 passing known planets 0.59–0.87, while 9 real signals at 1.5–2.3 passed, all
blends or EBs (TIC 304448456, found an EB on a neighbour by follow-up, among them).

**4. Secondary: ≥ 5σ near phase 0.5 (±0.1), ≥ 7σ elsewhere (was 5σ anywhere).** 5 bright
injections (expected SNR 20–355) were rejected for 5–6.7σ dips at phases 0.06–0.65; the
off-phase secondaries of sector-104 EBs are 9–15σ. New metric `secondary_half_snr`.

**5. Rotation on the raw flux: metric only, not a flag.** The variability test runs on the
flattened flux, where rotation slower than ~1–2 d is gone. A Lomb–Scargle of the raw flux (a
line removed per orbit segment, 0.1–8 d) now gives `rot_period` and `rotation_match`. As a soft
flag it fired on 20 % of all stars, on 18 of 342 passing injections and on TOI-177 b: raw TESS
flux is dominated by systematics at 4–8 d. It is shown on the report page instead.

**Effect on sector 104** (same SDE ≥ 10 / SNR_TLS ≥ 6.5): 31 signals pass (was 40); all 8
known planets still pass; the 9 removed are the duration cases above (6 "new", 3 known
neighbour EBs); none added. "New" candidates: 13 (was 19). Injections: 342 of 412 found pass
(was 339); completeness at SNR_exp ≥ 7: 79.9 % (was 79.2 %); losses to non-threshold tests:
16 (was 26).

**Follow-up, re-run with the fixes.** The five earlier verdicts stand (now with 6 of 7
repeating sectors for TIC 306101698; TIC 394585043 decided by the combined other sectors, since
no single one could show it). The six C-grade candidates never checked before: four are blends
on neighbours (424913827, 61615837, 88516648, 64336445), 372868160 does not repeat, and
214311908 repeats over six years, fits the star's density (0.84×), and stays ambiguous among
3–5″ Gaia neighbours (FPP 0.64). Details in `docs/REVIEW_sector104.md` and notebook 9 §7.

**Re-check.** `uv run pytest tests/test_vet.py tests/test_followup.py`; the comparison above
is in notebook 8 §6.1. Pre-fix tables are kept in `data/_before_fix/`.

## 2026-10-02 — Thresholds from injection-recovery: SDE ≥ 10, SNR_TLS ≥ 6.5

**Data.** 1,000 injections (`src/inject.py`, sector 104, P 0.3–12 d and Rp 0.5–6 R⊕
log-uniform, b uniform 0–0.9), re-vetted with the final `vet.py`. Analysis in
`notebooks/pipeline-8-injection-recovery.ipynb`.

**Pipeline sensitivity.** The search finds 50 % of injected planets at an expected SNR of 7.0
(plateau 94.6 %). After vetting, 50 % survive at an expected SNR of 9.5 (plateau 88.8 %). The plan
goal "> 80 % recovery at 3 mmag, P < 10 d" is **not met**: 72 % found and vetted, 80 % found
by the search alone.

**Threshold choice.** A grid of SNR column (TLS's / ours) × snr_min (5–8) × sde_min (8, 9, 10)
was scored by completeness (injections with expected SNR ≥ 7 that are found and pass) and by the
number of real sector-104 signals passing (almost all false positives, per the review). Two
settings dominate the textbook SDE ≥ 9 / SNR ≥ 7 (76.4 %, 42 real): SNR_TLS ≥ 6.5 with SDE ≥ 10
(**79.2 %, 40**, chosen as the most complete of the two) and our SNR ≥ 6.5 with SDE ≥ 10
(77.6 %, 39). Raising SDE from 9 to 10 costs injected planets almost nothing (≤ 1 %), because
injected planets that clear the SNR cut have high SDE anyway, while it removes real junk.

**Check against the by-eye review.** The new thresholds remove 4 of the 21 previously passing
"new" candidates, and 3 of the 4 were graded D (TIC 49143293, 308210382, 631919622; the
fourth, 389683926, was a C). They add one (TIC 64336445). All 8 known planets still pass
(lowest SDE 12.3, TOI-244 b). TOI-122 b still fails (SNR_TLS 5.6).

**Not chosen: our SNR.** The local-noise, red-noise-corrected SNR is more honest per signal
(TIC 49143293: 5.4 vs 8.5) but, at equal false-alarm cost, slightly less complete. It is kept as
a metric and selectable with `--snr-column snr`.

## 2026-10-02 — New hard test `transit_noise` (glitch-driven signals)

**Seen in the review** (`docs/REVIEW_sector104.md`, rank 7): TIC 394134586 passed all v2 tests
with a 1.3 R⊕ "planet". Its last three predicted transits each contain a ~1-h burst of strongly
increased scatter (+1 % to −2 %) whose mean is low. The seven earlier transits show nothing.

**Test.** For each transit, the robust scatter of the points in the central half of the
transit window, divided by the local out-of-transit scatter (±0.5 d). A transit shifts points
but keeps their scatter, so planets give ≈ 1. Transits above 2.0 are dropped and the signal must
survive (≥ 0.5 of its expected SNR, ≥ 2 transits left).

**Calibration** on sector 104: 8 known planets, max ratio 1.13–1.72 (none > 2); A/B candidates
≤ 1.72; TIC 394134586 has 3 transits > 2 (max 3.93). The central half avoids the ingress/egress
gradient, which would otherwise inflate the ratio for deep transits. Deep EBs (V-shaped) still
score high, which is harmless since other tests reject them. Unit test:
`test_noise_bursts_inside_transit_windows_fail`.

## 2026-10-02 — Vetting v2: four changes motivated by the by-eye review of the report

Reviewing the 32 "new" candidates page by page (`docs/REVIEW_sector104.md`) showed failure
modes the v1 tests missed. Each change is tested in `tests/test_vet.py`:

1. **Variability, harmonic matching in frequency.** A Lomb–Scargle peak is located only to
   ~1/baseline in frequency (~8 % in period at P = 2 d over 25 d). The v1 rule ("period ratio
   within 0.5 % of an integer") let obvious rotators through at 0.7 % and 1.3 % (TIC 32355300,
   TIC 183589961). Now: |f_var − n f_tls| (or |f_tls − n f_var|) < 1/baseline.
2. **Variability, significance.** v1 used any peak, however weak. TIC 161313862 (a clean-looking
   transit on a flat light curve) was rejected because of a power-0.003 peak (white-noise FAP
   ~2 %). Now the peak must have FAP < 10⁻³. New soft flag `variable_star`: a significant
   modulation ≥ 30 % of the depth at *any* period (TIC 160188205 has one at a non-harmonic
   period).
3. **New hard test `data_edges`.** Several "transits" sat at segment starts, right after gaps or
   just before momentum dumps (TIC 119624569: both transits 0.2–0.3 d before a dump, with data
   on one side only). Each transit now needs ≥ 50 % of the expected points on *both* sides
   (window max(duration, 1 h)), and the signal must survive dropping those that don't (same
   survival rule as the dump test, ≥ 2 transits left). In addition, the momentum-dump window
   grew from ±0.1 d to **±0.25 d**: the 0.5-d biweight smears a flux step over ±window/2, the
   same reasoning as `edge_trim`.
4. **Local noise per transit.** v1 gave every transit the global noise level. A transit in a
   stretch of strong scattered light (TIC 49143293) then looked as significant as one in quiet
   data. Now each transit's error uses the out-of-transit scatter within ±0.5 d (floored at half
   the global scatter). This SNR (`snr`, ours) is now often well below TLS's `snr_tls` for
   systematics-dominated signals (TIC 49143293: 5.4 vs 8.5). Which one the SNR threshold should
   use is decided with injection-recovery (notebook 8).

**Effect on sector 104:** 44 signals pass (was 52), and all 8 known planets still pass.
Injections are now re-run with stored light curves (`data/injection/lc/`), so future vetting
changes only need `python -m src.inject --revet` (minutes), not a new TLS run (hours).

## 2026-10-02 — Single-event test rejected every 2-transit signal (found by injection-recovery)

**Bug.** `_survival()` returned 0 whenever fewer than 2 transits were left after dropping the
suspect ones. For the single-event test, which drops the strongest transit, a 2-transit signal
always kept 1, so every planet with P ≳ 8 d failed. In the first 132 injections, 8 correctly
found planets (expected SNR 9–58) were rejected this way, and none were rejected for any other
reason.

**Fix.** `min_keep = 1` for the single-event test (the weaker transit must carry half its
expected SNR) and `min_keep = 2` for the momentum-dump test (two transits must lie away from the
dumps). Regression test: `test_two_transit_planet_passes_single_event`.

**Effect on sector 104.** 52 signals pass vetting (was 49). The 18 injections with 2 transits
were re-run with the fixed code.

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
