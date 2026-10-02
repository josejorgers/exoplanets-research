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
| 1 | 394585043 | 4.496 | 7,335 | 2.4 | 6.6 | – | **D** | **Not repeating**: searched in 9 other sectors (2018–2026, 51k points) within ±2 % of P, SDE ≈ 6 (noise). A 7,300-ppm transit would be obvious there. Notebook 9 §2. |
| 2 | 49143293 | 9.049 | 42,592 | 6.0 | 5.4 | – | **D** | Systematics: QLP flux jumps 0.4→2.0; the three "transits" sit at a segment start, on the 4185.04 dump, and in the scattered-light stretch after the gap. Local-noise SNR 5.4 vs TLS 8.5. |
| 3 | 229804549 | 0.732 | 2,385 | 2.8 | 17.3 | centroid | **B** | 30 consistent transits, clean periodogram with harmonics, flat odd/even. Concerns: dip wings broader than the model, a significant modulation at 14 % of the depth (just below the 15 % rejection limit), centroid 6σ. |
| 4 | 306101698 | 2.879 | 4,882 | 3.1 | 10.5 | centroid | **A** | Quiet star, box-shaped 2-h transit, 8 transits, odd/even 0.06σ, no secondary, clean periodogram peak with aliases. Only concern: centroid shift 8.6σ (possible blend; check pixels). |
| 5 | 61615837 | 1.108 | 3,558 | 2.4 | 10.3 | centroid | C | Large ramps in the raw flux; shallow, broad fold. |
| 6 | 304448456 | 0.592 | 4,357 | 2.0 | 9.7 | centroid | **B** | Box-shaped dip, SDE 15, clean periodogram; ramps are flattened well. Centroid 7.6σ, contratio 0.29. |
| 7 | 394134586 | 2.212 | 1,341 | 1.3 | 9.5 | contamination | **D** | Glitch-driven: three 1–2-cadence downward spikes (−1 to −2 %) in orbit 2 coincide with the predicted times; much shorter than the 1.8-h duration; nothing in orbit 1. |
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

**Summary:** 1 A, 2 B, 11 C, 9 D. The A and B candidates (TIC 306101698, 229804549,
304448456), plus TIC 161313862 (whose centroid issue should be settled one way or the other)
and TIC 394585043 (as a worked example of a candidate that does not repeat), go to notebook 9.

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

- **Glitch-driven signals** (ranks 7, 13): a few 1–2-cadence spikes much shorter than the
  transit duration. A test comparing the in-transit point distribution with a transit shape
  (e.g. the fraction of the depth carried by the lowest 10 % of in-transit points) would
  catch them.
- **Strong ramps** (C grades): QLP light curves with 10–50 % trends still yield dips after
  flattening. Comparing per-transit depths with the local slope of the raw flux could quantify
  this.
- **Centroid shifts** are frequent among plausible candidates (6–16σ). Only pixel-level
  difference imaging can say *where* the signal comes from (notebook 9 §4).
