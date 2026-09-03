# Annex A. The twenty-nine compared quantities behind Figure 11

*Prepared 2026-09-02 for the revision of "Boundary Condensation of a Two-Layer Dirichlet-to-Neumann Operator: Trace Conformity and Air-Gap Mesh Reduction".*

Figure 11 of the submitted manuscript reports two counts a reader cannot check from the plot alone: *twenty-one of twenty-nine fall within 5 %*, and *twenty of twenty-six* once the three quantities the manuscript declares not validated are set aside. This annex tabulates the twenty-nine bars at full precision, with the model value, the reference value, the signed deviation and the transcript each is read from, so that a referee can recount both figures line by line.

**Convention.** The deviation is the signed relative deviation of (38), formed on full-precision values and rounded once for display; negative means the network is below the reference. The threshold is strict, `|deviation| < 5`. Rows 25 and 29 are efficiencies: Section 5.3 requires them in percentage points, and the two columns for them are given below.

**Provenance.** Column *source* names the transcript each line is read from:

- `R7` = `outputs/MEC_BLDC/R7_scorecard_out.txt`, produced by `code/MEC_BLDC/RUN_R7_SCORECARD.m`;
- `G4` = `outputs/MEC_BLDC/G4_points_out.txt`, produced by `code/MEC_BLDC/RUN_G4_POINTS.m`, which re-reads `code/MEC_BLDC/R7_scorecard.mat`.

Every value below was extracted mechanically from those two files by `code/revision/build_manifest.py`; none was transcribed by hand.

## A.1 The twenty-nine quantities

| # | quantity | model | reference | deviation (%) | within 5 % | source |
|---:|---|---:|---:|---:|:---:|:---|
| 1 | Fundamental air-gap flux density B_g1 (T) | 1.07755 | 1.07455 | +0.2788 | yes | R7 |
| 2 | Mean air-gap flux density |B_r| (T) | 0.761473 | 0.760909 | +0.0741 | yes | R7 |
| 3 | Peak air-gap flux density B_r (T) | 0.984576 | 0.990727 | -0.6208 | yes | R7 |
| 4 | Tangential air-gap flux density B_t, rms (T) | 0.143165 | 0.133238 | +7.4505 | **no** | R7 |
| 5 | Flux-line value A, peak (Wb/m) | 0.00585972 | 0.00585343 | +0.1074 | yes | R7 |
| 6 | Magnetomotive force per magnet, mean (A) | 763.243 | 769.445 | -0.8061 | yes | R7 |
| 7 | Magnetomotive force per gap, mean (A) | 723.829 | 721.504 | +0.3223 | yes | R7 |
| 8 | Magnet magnetomotive force, dispersion (%) | 6.8073 | 6.58458 | +3.3826 | yes | R7 |
| 9 | Reluctance factor k_r (-) | 1.05445 | 1.08693 | -2.9880 | yes | R7 |
| 10 | Leakage factor k_l (-) | 0.853766 | 0.86723 | -1.5525 | yes | R7 |
| 11 | Self inductance L_a (mH) | 50.8159 | 50.2086 | +1.2096 | yes | R7 |
| 12 | Mutual inductance M (mH) | -2.22152 | -2.13292 | +4.1538 | yes | R7 |
| 13 | Synchronous inductance L_d (mH) | 53.0374 | 52.3415 | +1.3295 | yes | R7 |
| 14 | Phase back-electromotive force, peak (V) | 240.936 | 233.031 | +3.3925 | yes | R7 |
| 15 | Six-step envelope, peak (V) | 469.262 | 459.918 | +2.0315 | yes | R7 |
| 16 | Flux linkage, peak (Wb) | 0.261206 | 0.258328 | +1.1139 | yes | R7 |
| 17 | Electromagnetic torque on load (N m) | 4.94901 | 4.86664 | +1.6925 | yes | R7 |
| 18 | Phase current, rms, on load (A) | 1.41857 | 1.55162 | -8.5752 | **no** | R7 |
| 19 | Iron loss at no load (W) | 21.3358 | 19.7298 | +8.1401 | **no** | R7 |
| 20 | Magnet eddy-current loss at no load (W) | 0.212037 | 0.333466 | -36.4143 | **no** | R7 |
| 21 | Settled speed on load (rpm) | 1257.03 | 1339.83 | -6.1797 | **no** | R7 |
| 22 | Mean dc-bus current (A) | 1.42741 | 1.65746 | -13.8796 | **no** | R7 |
| 23 | Iron loss on load (W) | 19.083 | 18.465 | +3.3471 | yes | R7 |
| 24 | Magnet eddy-current loss on load (W) | 1.57752 | 3.18268 | -50.4342 | **no** | R7 |
| 25 | Efficiency on load (%) | 86.334 | 86.5148 | -0.2089 | yes | R7 |
| 26 | Back-electromotive-force constant k_E (mV/rpm) | 312.841 | 323.644 | -3.3380 | yes | R7 |
| 27 | Standstill torque, speed sweep (N m) | 30.6807 | 29.2331 | +4.9520 | yes | R7 |
| 28 | Peak power, speed sweep (W) | 1115.07 | 1495.27 | -25.4270 | **no** | R7 |
| 29 | Peak efficiency, speed sweep (%) | 85.4456 | 89.4319 | -4.4574 | yes | R7 |

## A.2 The two efficiency rows in percentage points

Section 5.3 reports a deviation between two efficiencies in points, not in per cent. `G4` carries both readings at full precision.

| # | quantity | model (full) | reference (full) | relative (%) | points (p.p.) | source |
|---:|---|---:|---:|---:|---:|:---|
| 25 | Efficiency on load (%) | 86.3340330680 | 86.5147786624 | -0.208919 | -0.180746 | G4 |
| 29 | Peak efficiency, speed sweep (%) | 85.4455763512 | 89.4319268475 | -4.457413 | -3.986350 | G4 |

## A.3 The three quantities declared not validated

The manuscript excludes three quantities from the second count and says why. The motives below are the transcript's own.

| # | quantity | deviation (%) | motive as recorded in R7 |
|---:|---|---:|:---|
| 20 | Magnet eddy-current loss at no load (W) | -36.4143 | §6.6 : "The reference cannot resolve the quantity on which it is being compared." |
| 24 | Magnet eddy-current loss on load (W) | -50.4342 | §6.6 + Table 19 : ecart -1,6 W sous le plancher de resolution de 52,7 W |
| 27 | Standstill torque, speed sweep (N m) | +4.9520 | §5 : "cannot be claimed as a validation ... excluded from the comparison" |

## A.4 Re-derivation of the two published counts

Counted directly from the table of A.1, with the strict threshold `|dev| < 5`:

| count | re-derived here | printed by R7 | printed by G4 | published in Fig. 11 | holds |
|:--|:--|:--|:--|:--|:--|
| within 5 %, all twenty-nine | **21 / 29** | 21 / 29 | 21 / 29 | twenty-one of twenty-nine | yes |
| within 5 %, excluding the three not validated | **20 / 26** | 20 / 26 | 20 / 26 | twenty of twenty-six | yes |

The quantities outside the threshold are, in order: #4 Tangential air-gap flux density B_t, rms (T) (+7.4505 %), #18 Phase current, rms, on load (A) (-8.5752 %), #19 Iron loss at no load (W) (+8.1401 %), #20 Magnet eddy-current loss at no load (W) (-36.4143 %), #21 Settled speed on load (rpm) (-6.1797 %), #22 Mean dc-bus current (A) (-13.8796 %), #24 Magnet eddy-current loss on load (W) (-50.4342 %), #28 Peak power, speed sweep (W) (-25.4270 %).

Removing #20, #24, #27 from that list leaves the 6 exceedances of the second count.

**Quantities within one point of the threshold.** In the relative convention there are 3: #12 (+4.1538 %), #27 (+4.9520 %), #29 (-4.4574 %). G4 records "3 en relatif, 2 en points". The submitted caption says *two*, which is the reading in the points convention of Section 5.3, where row 29 becomes -3.986350 p.p. and leaves the one-point band; the caption and the transcript agree once the convention is stated.

**Guard.** G4 recomposes the twenty-nine deviations from the saved model and reference values and compares them with the stored deviation column; the maximum discrepancy is 7.105e-15, that is machine precision.

## A.5 What this annex does not settle

Two limits belong with the table.

1. The reference carries no uncertainty band on any field quantity, and the manuscript states that it does not. A deviation of a few tenths of a per cent on rows 1 to 8 is therefore not evidence of agreement to that accuracy.
2. Rows 17 to 25 are read at an operating point the reference does not share: the network settles at a different speed. The efficiency agreement of row 25 is obtained there and nowhere else.

---

Files read to build this annex: `outputs/MEC_BLDC/R7_scorecard_out.txt`, `outputs/MEC_BLDC/G4_points_out.txt`.
