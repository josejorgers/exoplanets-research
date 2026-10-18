# By-eye review of the sector-104 candidates

Report reviewed: `data/reports/sector_104.pdf` (vetting v2, 2026-10-02), with pages also in
`data/reports/sector_104/{rank}_{tic}.png`. How to read a page:
`notebooks/pipeline-7-report-walkthrough.ipynb`.

## Method

Every signal that passes all hard vetting tests and is not explained by a known object
(`known_class = new`) was looked at page by page. The questions, in order (notebook 7, §5):

1. **Raw light curve (panel 1).** Is the star quiet, variable, or full of systematics (ramps,
   steps, glitches)? Do the transits fall where the raw flux misbehaves?
2. **Flattened light curve (panel 2).** Are the individual transits visible or at least not
   contradicted? Do they avoid gaps, segment starts and momentum dumps?
3. **Fold (panel 3).** Is there a transit-shaped dip that the binned points follow, with flat
   wings? Or is it a few low points, a broad V, or a step?
4. **Odd/even and phase 0.5 (panels 4, 5).** Consistent depths? Nothing at phase 0.5?
5. **Periodogram (panel 6).** A sharp peak with aliases, or one peak in a forest?
6. **Table.** Soft flags (centroid, contamination, V shape, variable star) and how strong they
   are.

**Verdict scale**

| grade | meaning | next step |
|---|---|---|
| **A** | looks like a planet in every panel; at most one soft flag | full follow-up (notebook 9) |
| **B** | plausible; one clear concern | follow-up |
| **C** | weak or ambiguous: low SNR, systematics nearby, several flags | keep; revisit with more sectors |
| **D** | not a planet by eye: systematics, variability, glitches, likely blend | reject (reason recorded) |

## Candidates (`known_class = new`, passes all hard tests)

| rank | TIC | P [d] | depth [ppm] | Rp [R⊕] | SNR | flags | grade | reason |
|---|---|---|---|---|---|---|---|---|
| 1 | 394585043 | 4.496 | 7,335 | 2.4 | 6.6 | – | **D** | **Not repeating**: the 9 other sectors (2018–2026, 51k points) searched together within ±2 % of P return a different period at a third of the depth, where this transit would reach SNR ≈ 16. (No single other sector is sensitive enough on its own.) Notebook 9 §1. |
| 2 | 49143293 | 9.049 | 42,592 | 6.0 | 5.4 | – | **D** | Systematics: QLP flux jumps 0.4→2.0; the three "transits" sit at a segment start, on the 4185.04 dump, and in the scattered-light stretch after the gap. Local-noise SNR 5.4 vs TLS 8.5. |
| 3 | 229804549 | 0.732 | 2,385 | 2.8 | 17.3 | centroid | **B** | 30 consistent transits, clean periodogram with harmonics, flat odd/even. Concerns: dip wings broader than the model, a significant modulation at 14 % of the depth (just below the 15 % rejection limit), centroid 6σ. |
| 4 | 306101698 | 2.879 | 4,882 | 3.1 | 10.5 | centroid | **A** | Quiet star, box-shaped 2-h transit, 8 transits, odd/even 0.06σ, no secondary, clean periodogram peak with aliases. Only concern: centroid shift 8.6σ (possible blend; check pixels). |
| 5 | 61615837 | 1.108 | 3,558 | 2.4 | 10.3 | centroid | C | Large ramps in the raw flux; shallow, broad fold. |
| 6 | 304448456 | 0.592 | 4,357 | 2.0 | 9.7 | centroid | **B** | Box-shaped dip, SDE 15, clean periodogram; ramps are flattened well. Centroid 7.6σ, contratio 0.29. |
| 7 | 394134586 | 2.212 | 1,341 | 1.3 | 9.5 | contamination | **D** | Glitch-driven: the last three predicted transits each contain a ~1-h burst of strongly increased scatter (+1 % to −2 %) whose mean is low; the seven earlier transits show nothing. |
| 8 | 308210382 | 10.947 | 15,565 | 3.2 | 8.0 | contamination | **D** | Two narrow events, one next to the 4198.2 dump; raw flux shows glitch-like drops. |
| 9 | 161313862 | 0.293 | 2,705 | 1.2 | 7.9 | centroid | C | Shallow 7-h-period signal; 4.3σ dip at phase 0.52 (threshold 5); centroid shift 16σ → probably a blend. |
| 10 | 389683926 | 9.654 | 9,881 | 3.6 | 7.3 | contamination | C | Three transits on strong ramps; contratio 0.65. |
| 11 | 88516648 | 1.214 | 847 | 1.8 | 7.2 | centroid | C | Near the SNR limit; fold barely distinguishable from noise. |
| 12 | 424913827 | 1.573 | 676 | 0.9 | 7.1 | centroid | C | Sub-Earth-size signal at the SNR limit on a star with ramps. |
| 13 | 231888176 | 5.314 | 2,943 | 1.7 | 6.4 | contamination | **D** | Short glitch-like drops in the raw flux of orbit 2 coincide with transits (same pattern as rank 7). |
| 14 | 160188205 | 0.701 | 8,511 | 4.2 | 14.2 | centroid, variable_star | **D** | Strongly variable star; broad V-shaped "dip" that is the residual of the modulation (dominant period 0.13 d, not a harmonic of P). |
| 15 | 218663210 | 0.788 | 7,856 | 2.9 | 12.1 | v_shaped, centroid | **D** | V-shaped (T23/T14 = 0.10), centroid shift 59σ, contratio 0.49: most likely a blended eclipsing binary. |
| 16 | 217484923 | 9.006 | 2,789 | 1.5 | 10.6 | centroid, contamination | C | Three long (3.9 h) transits on ramps. |
| 17 | 214311908 | 5.692 | 13,793 | 4.4 | 10.2 | centroid, contamination | C | Three transits; very strong ramps; contratio 2.2. |
| 18 | 372868160 | 4.513 | 2,053 | 1.5 | 7.6 | centroid, contamination | C | Weak; contratio 1.0. |
| 19 | 213062150 | 1.685 | 4,902 | 1.5 | 7.5 | v_shaped, centroid | C | V-shaped, weak. |
| 20 | 74510850 | 3.048 | 3,831 | 1.5 | 7.5 | centroid, contamination | C | Ramps; contratio 2.0. |
| 21 | 631919622 | 7.565 | 11,192 | 5.3 | 5.9 | v_shaped, contamination | **D** | Large irregular variability (±5 %); the dip is a residual of it. |
| 22 | 389127976 | 0.576 | 30,288 | 5.9 | 7.1 | centroid, contamination, variable_star | **D** | ~20 % flux steps in the raw data; contratio 7.2; the fold is a step, not a transit. |
| 23 | 470329609 | 0.854 | 1,019 | 1.0 | 7.5 | 4 flags | C | Weak, V-shaped, every flag raised. |

**Summary:** 1 A, 2 B, 11 C, 9 D (23 signals: the "new" list of vetting v2, before the
`transit_noise` test and the threshold tuning below). The A and B candidates (TIC 306101698, 229804549,
304448456), plus TIC 161313862 (whose centroid issue should be settled one way or the other)
and TIC 394585043 (as a worked example of a candidate that does not repeat), went to notebook 9.
None survived (see "Follow-up" below).

## After threshold tuning (injection-recovery, SDE ≥ 10 and SNR ≥ 6.5)

The `transit_noise` test (added after this review) removes rank 7, TIC 394134586, and the
thresholds chosen in notebook 8 change the list further (19 candidates instead of 23):

| TIC | change | grade | comment |
|---|---|---|---|
| 49143293 | removed (SDE 9.95) | D | systematics, see above |
| 308210382 | removed (SDE 9.71; also fails `transit_noise`) | D | glitch-like events |
| 631919622 | removed (SDE 9.31) | D | irregular variability |
| 389683926 | removed (SDE 9.09) | C | three transits on ramps |
| 64336445 | **added** (SNR_TLS 6.96, below the old 7.0) | C | noisy light curve (~1.5 % scatter) with strong ramps; centroid 9σ, contratio 0.56. Probably a blend or noise. |

Three of the four removed candidates had been graded D by eye. That is consistent with the
tuned thresholds removing junk rather than planets, but it is weak evidence: the thresholds
were chosen partly to minimise the number of real sector-104 signals passing, almost all of
which are junk, and with 9 of 23 graded D, three D's among four removals happens by chance
about one time in six.

## Second-pass search (a second planet around the candidates)

`search --second-pass --only-candidates` masked each of the 40 passing signals and searched
again. 7 stars have a second signal with SDE > 9, and none passes vetting:

- TIC 229804549 returns at the **same** period (masking 2× the TLS duration does not cover the
  dip's broad wings). This supports the "broader than a transit" concern.
- TIC 262530407 (GJ 3090, host of TOI-177 b): SDE 15.4 at 6.56 d, failing `momentum_dump`,
  `data_edges` and `secondary`. These are systematics, and not GJ 3090 c (P = 15.9 d, known
  from radial velocities, outside a one-sector search).
- The rest are rejected by SDE/SNR, or are aliases of EB signals on known-neighbour stars.

## Follow-up of the five best candidates (notebook 9)

Every sector of each star was searched again (narrow window around P), the pixels were
difference-imaged, the combined light curve was fitted with an MCMC transit model, and
TRICERATOPS computed a false-positive probability. LHS 3844 b went through the same chain as
a control.

| TIC | grade | repeats in other sectors | pixels | fitted ρ★ / TIC ρ★ | FPP | outcome |
|---|---|---|---|---|---|---|
| LHS 3844 (control) | known planet | 10 of 10 | on target (0.04 ± 0.08 px) | 0.80 | 0.015 | **planet**; FPP at the validation limit (optimistic: background scenarios were skipped, see below) |
| 306101698 | A | 6 of 7 | **off target**: loss 2 px away | 0.07 | 0.90 (NEB) | **blended EB on a neighbour** |
| 304448456 | B | 2 of 2 at P but 2× deeper; deeper again at 2P | **off target**: loss 2 px away | 0.04 | 0.99 | **EB on a neighbour, found at half its period** |
| 229804549 | B | 2 of 7 (depth 45 %) | **off target**: loss 6 px away | 0.04 | 1.00 | **blend; does not repeat cleanly** |
| 161313862 | C | 0 of 4 informative (6 observed) | inconclusive | – | 0.27 | **does not repeat** |
| 394585043 | D | no single sector informative; 9 combined: absent | inconclusive | – | 0.06 | **does not repeat** (was the pipeline's #1) |

(Repetition column after the 2026-10-04 fix: the first version's "narrow-window" searches were
full searches in most sectors, and counted sectors too weak to show the signal. The verdicts
did not change; see `docs/NOTES_decisions.md`.)

All FPPs in this table were computed without TRICERATOPS' background-star scenarios (the
TRILEGAL server was unreachable), so they are lower limits. That does not change any outcome:
every candidate is rejected by repetition or pixels, independently of FPP.

**Result: no new planet candidate from sector 104 survives follow-up; no CTOI submission.**
The dominant false positive among vetting survivors is a neighbouring eclipsing binary, which
vetting on the target's own light curve cannot exclude. The repetition check and the difference
image are cheap enough to run on every survivor, so they are the natural next additions to the
automated pipeline.

## Signals rejected by vetting but looked at anyway

| TIC | rejected by | by eye | comment |
|---|---|---|---|
| 231702397 (TOI-122 b) | snr (TLS 5.6) | real planet | The SNR threshold costs this known planet; see notebook 8 for the trade-off. |
| 119624569 | momentum_dump, data_edges | agree | Both "transits" 0.2–0.3 d before dumps, with data on one side only. Passed v1. |
| 397113523 | data_edges | agree | "Transits" at segment starts. Passed v1. |
| 32355300, 183589961 | variability | agree | Obvious rotators (P_var ≈ P_tls); passed v1 because of a too-tight period-ratio tolerance. |
| 139722169 | variability | agree, weakly | Broad, sinusoid-like dip at a significant harmonic. Could also be a grazing EB; not a convincing planet either way. |
| 87329149 | variability | agree | Noisy; the dip is not transit-shaped. |

## What the review changed in the pipeline

The first review (vetting v1, 32 "new" candidates) found recurring failure modes that led to
vetting v2: the data-edges test, the ±0.25-d dump window, local per-transit noise, and
frequency-based variability matching with a significance requirement. See
`docs/NOTES_decisions.md` (2026-10-02) and notebook 5 §12.1.

Still open, seen in this review:

- **Glitch-driven signals** (ranks 7, 13): short bursts of increased scatter inside a few
  transit windows, whose mean happens to be low. A test comparing the in-transit scatter with
  the local out-of-transit scatter (a transit shifts the points down but keeps their scatter)
  would catch them.
- **Strong ramps** (C grades): QLP light curves with 10–50 % trends still yield dips after
  flattening. Comparing per-transit depths with the local slope of the raw flux could quantify
  this.
- **Centroid shifts** are frequent among plausible candidates (6–16σ). Only pixel-level
  difference imaging can say *where* the signal comes from (notebook 9 §4).

## Vetting v3 and follow-up of every remaining candidate (2026-10-04)

Vetting v3 (`docs/NOTES_decisions.md`, 2026-10-04) tightened the duration test (≤ 1.5× the
central-transit duration, no lower bound) and made the secondary test phase-aware. It removes
six "new" signals, all graded C or worse here or rejected by follow-up: TIC 304448456 (B, EB
on a neighbour), 160188205 (D), 217484923, 213062150, 74510850, 470329609 (C). Nothing was
added. **13 "new" candidates remain.**

The six remaining C grades, never checked before, went through the repetition check and the
difference image; the two that repeated went through the full chain (fit, FPP).

| TIC | P [d] | repeats in other sectors | pixels | outcome |
|---|---|---|---|---|
| 424913827 | 1.573 | 2 of 4, depth up to 3.6× larger (sector 94: SNR 26); deeper at 2P | **off target**: S/N 200 loss 3.6 px away | **EB on a neighbour, found at half its period** |
| 61615837 | 1.108 | 1 of 3 informative, 3× deeper (sector 67: SNR 24) | **off target**: loss 2.7 px away | **EB on a neighbour** |
| 88516648 | 1.214 | 1 of 3 informative | **off target**: loss 2.4 px away | **blend** |
| 64336445 | 0.677 | only in all sectors combined (depth 0.6×) | inconclusive (S/N < 5); light centre moves 2.4 px | **likely blend**: fitted stellar density 0.10× the target's (FPP 0.45) |
| 372868160 | 4.513 | 0 of 4 informative; combined absent | – | **does not repeat** |
| 214311908 | 5.692 | **3 of 5** (2 of 2 informative), same depth; combined too | loss 0.3 px from the catalog position | **ambiguous**, see below |

**TIC 214311908** passes everything the data can test: it repeats over six years (2019–2026)
at the same 13,800-ppm depth, the MCMC fit gives Rp = 4.9 ± 0.4 R⊕, b ≈ 0.6 and a stellar
density 0.84× the TIC value, the transit is U-shaped (T23/T14 = 0.58), odd/even agree (0.9σ),
and nothing is seen at phase 0.5. But the field is crowded (contratio 2.2; ~11 % of the aperture
light is the target's), and Gaia DR3 lists stars at 3.2″ (ΔG = 3.2), 5.2″ and 5.4″ that are
bright enough to make the dip if they were eclipsing binaries, unresolvable at 21″ per pixel.
TRICERATOPS gives FPP = NFPP = 0.64 (top scenario: planet on the target, 33 %), and that is
optimistic (background scenarios skipped). It does not meet our FPP < 0.5 bar for a CTOI.
Deciding it needs seeing-limited photometry of the 3–5″ neighbours during a predicted transit.

(Its difference image was first labelled "off target": the flux loss is 0.76 px from the
out-of-transit light centroid, which the neighbours pull away from the star. The off-target
rule now requires the loss to be away from the catalog position as well; no other verdict
changed.)

**Updated result for sector 104:** 12 of the 13 "new" candidates are explained (7 eclipsing
binaries or blends on neighbours counting the V-shaped TIC 218663210, 3 that do not repeat, 2
systematics); one, TIC 214311908, is an unresolved candidate that TESS data alone cannot
settle. No CTOI submission.

