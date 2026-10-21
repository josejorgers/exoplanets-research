# CTOI draft: TIC 214311908

**Status: draft, not submitted.** ExoFOP accepts community candidates only after they are published
in a peer-reviewed journal, or in an AAS Research Note (RNAAS) that cites a detection and vetting
method already published in a peer-reviewed journal. The paper URL is a required field
(https://exofop.ipac.caltech.edu/tess/candidate_help.php, checked 2026-10-21). So this file is
material for that paper, or for a collaborator's, not something to upload now. Re-check the target page (https://exofop.ipac.caltech.edu/tess/target.php?id=214311908) on the day you submit:
on 2026-10-19 it listed no TOI, no CTOI, no imaging, spectroscopy, time series or files.

Numbers come from `ctoi_parameters.json` (written by `scripts/ctoi_TIC214311908.py` from the
MCMC fit of all six sectors and the TRICERATOPS runs in `data/followup/214311908/`).

## Form fields

| field | value |
|---|---|
| TIC ID | 214311908 |
| Period (days) | 5.691498 ± 0.000007 |
| Epoch (BJD_TDB) | 2460852.2989 ± 0.0012 |
| Depth (ppm) | 15100 ± 800 (fitted mid-transit depth; QLP light curves are deblended) |
| Duration (hours) | 1.71 ± 0.10 (first to fourth contact) |
| Planet radius (R⊕) | 4.85 ± 0.42 (TIC R★ = 0.381 R☉) |
| Insolation (S⊕) | 7.7 ± 2.8 |
| Equilibrium temperature (K) | 464 ± 39 (albedo 0, full redistribution) |
| SNR | 19 (expected, all six sectors; discovery sector alone 9.9) |
| Proposed disposition | PC |
| Sectors used | 12, 39, 66, 93, 103, 104 (QLP FFI light curves; no 2-min data exist) |
| User / group / tag | yours to fill |

The epoch is the fit's reference epoch, near the middle of the data. The 1σ uncertainty on
predicted transit times through April 2027 is about 2.3 minutes
(`data/followup/214311908/ephemeris.csv`).

## Notes field (short version)

> Found in an independent TLS search of nearby M dwarfs in TESS S104 (QLP), recovered in 5 other
> sectors (2019-2026) at consistent depth. U-shaped; odd/even consistent; no secondary; fitted
> stellar density 0.94x TIC. Crowded field (TIC contratio 2.2). PRF fit of the difference images
> (S39-S104) places the source within 0.4 px of the target and excludes all resolved Gaia
> neighbours >6" (incl. a co-moving companion at 11.7" and a Mira at 17.8"). The 3.1" (dG=3.2) and
> 4.3"/5.7" stars are disfavoured by the transit shape (3.1": ~3.5 sigma). The two Gaia DR3 EBs nearby (26.7",
> P=0.316 d; 62", P=0.261 d) are too faint and at unrelated periods. TRICERATOPS (incl. TRILEGAL)
> FPP 0.09 +/- 0.05 with pixel-excluded stars removed; 0.45 +/- 0.14 without. No high-resolution
> imaging yet. Suggested follow-up: seeing-limited photometry with a <=2" aperture (an EB on the
> 3.1" star would show a ~40% eclipse) and high-resolution imaging.

ExoFOP may limit the length of this field: check the form, and if it is too long put the short
version in the field and upload this file as a supporting document.

## Full description (for a supporting file or for TFOP)

**Detection.** TIC 214311908 (Tmag 13.01, M dwarf: Teff 3301 K, R★ 0.381 R☉, M★ 0.369 M☉,
d = 64.2 pc) was found by our pipeline in sector 104 (QLP 200-s FFI light curve, wotan biweight
0.5 d, TLS; SDE 14.3, SNR 9.9). It passed automatic vetting (odd/even, secondary, duration
against the stellar density, data edges, momentum dumps, transit-to-transit noise).

**Repetition.** Searched sector by sector (TLS restricted to ±2 % of P), it recurs at the same
period in sectors 39, 93 and 103. Sectors 93 and 103 are the only other sectors where the transit
should reach an expected SNR of 7.5, and both show it. Sector 39 also shows it, at expected SNR
7.1. Sectors 12 (30-min cadence) and 66 are too weak alone. A search of all sectors together
recovers P with SDE 24 (expected SNR ≈ 19), at a depth 0.95 times the discovery sector's
(figure 1). A search at twice the period finds no depth difference between alternate transits.

**Transit fit.** batman + emcee fit of all six sectors: Rp/R★ = 0.117 ± 0.006, b = 0.54
(+0.20/−0.36), a/R★ = 25 ± 4, T14 = 1.71 ± 0.10 h, Rp = 4.85 ± 0.42 R⊕. The stellar density from
the transit is 0.94 times the TIC value. The flat bottom limits the host to at most 3.1 mag
fainter than the target (99.7 %, Vanderburg et al. 2019 test).

**Neighbours and pixels (figures 2, 3).** The field is 8.5° from the Galactic plane, and the TIC
contamination ratio is 2.2. Every Gaia DR3 source within 60″ was considered:

- **Too faint:** most sources would need an eclipse deeper than 100 % to make the dip.
- **Excluded by the PRF fit:** we fitted the in-minus-out-of-transit difference images of
  sectors 39–104 with the TESS PRF. Every neighbour bright enough and farther than 6″ is worse
  than the target as the source at Δχ² > 25. This includes:
  - the co-moving companion at 11.7″ (G 14.65; RUWE 13.7, so itself probably a close binary),
    at Δχ² ≈ 86;
  - a long-period variable / Mira at 17.8″, at Δχ² ≈ 209;
  - a G = 17.0 star at 12.0″, at Δχ² ≈ 113.
- **Free-position fit:** the best-fit position of the lost flux is 0.07–0.37 px from the target
  in every sector.
- **The three faint stars within 6″:**
  - The 3.1″ star (G 17.5) can't be separated from the target by the pixels (Δχ² ≈ 5). It would
    need a ~30–50 % eclipse (G or T band), and the transit shape disfavours it at about 3.5σ.
  - The 4.3″ star (G 18.6) would need an ~80 % eclipse, and the shape excludes it.
  - The 5.7″ star (G 17.7) is excluded by both tests (Δχ² ≈ 27, and the shape).
- **Gaia data on the target:** it looks single (RUWE 1.08, no multi-peak flags, no
  non-single-star solution).
- **Known eclipsing binaries nearby (Gaia DR3, listed by SIMBAD):** neither can produce the
  signal. They are too faint to make a 1.5 % dip even when totally eclipsed, and their eclipses
  would recur many times per 5.69 d.
  - G = 19.50 at 26.7″, P = 0.3161 d;
  - G = 20.00 at 62.1″, P = 0.2609 d.

  Two more Gaia EBs, at 167″ and 168″ (P = 0.750 d and 1.595 d), are 8 px away. That is outside
  the aperture, where the pixel fit rejects every source, and their periods are unrelated.
- **Archival catalogues:** 2MASS and SkyMapper DR4 show no relevant source that Gaia lacks.

**False-positive probability.** TRICERATOPS was run five times, with:
- TRILEGAL background populations;
- Gaia-based parameters for neighbours that lack TIC values;
- TIC duplicates removed.

| run setting | FPP | NFPP |
|---|---|---|
| neighbours excluded by the PRF fit removed (4 TIC entries) | 0.09 ± 0.05 | 0.05 ± 0.03 |
| no neighbours removed | 0.45 ± 0.14 | |

**Open scenarios:**
- an eclipsing binary on the 3.1″ star;
- an unresolved companion closer than about 1″, which needs high-resolution imaging.

Seeing-limited photometry with an aperture of 2″ or less would test the first scenario. It should
show a 1.5 % transit on the target, or a roughly 40 % eclipse on the 3.1″ star. Five transits are
fully observable from CTIO, SAAO or Siding Spring between February and April 2027: 19 Feb, 25 Mar,
5 Apr, 11 Apr and 17 Apr (`data/followup/214311908/ephemeris.csv`).

**Products.**
- `fig1_lightcurve.png`
- `fig2_pixels.png`
- `fig3_neighbours.png`
- `summary.pdf`: the three figures together
- `ctoi_parameters.json`
- Analysis notebook: `notebooks/candidate-TIC214311908.ipynb` in
  https://github.com/josejorgers/exoplanets-research
