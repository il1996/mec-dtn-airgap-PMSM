# Independent numerical verification

**Manuscript** — *Boundary Condensation of a Two-Layer Dirichlet-to-Neumann Operator: Trace
Conformity and Air-Gap Mesh Reduction*, submitted to *Engineering Analysis with Boundary Elements*.
**Date** — 3 September 2026. **Scope** — every published value that the manuscript's own equations
determine.

## What this is, and what it is not

This is a second implementation, in Python, of the equations printed in the manuscript. It shares
no code path with the MATLAB chain under `code/`, reads nothing from `outputs/`, and needs neither
MATLAB nor a finite-element licence. It exists so that a reader can check the closed forms and the
tables that follow from them without the original chain.

It is **not** a reproduction of the coupled model. Anything that requires the reluctance network,
the Newton solve or the finite-element reference is outside it, and is marked as such below.

```
cd validation/independent
python3 verify_all.py            # the report
python3 verify_all.py --csv      # and results/independent_verification.csv
```

Requires Python 3.8 or later. The standard library only; no NumPy, no SciPy. Runtime ≈ 13 s.

| file | what it is |
|---|---|
| `dtn_operator.py` | the two-layer operator, rebuilt from (10), (13), (14), (26), (29)–(33) and the geometry of Table 3 |
| `published.py` | a transcription of the printed page — the values being checked against |
| `verify_all.py` | the checks and the report |
| `results/verify_all_out.txt` | the transcript of the run |
| `results/independent_verification.csv` | one row per check, machine-readable |

## Result

```
checks: 251   agreeing: 248   mismatches: 0
not reproducible from the manuscript alone: 1   skipped for cost: 2
VERDICT: PASS
```

| verdict | count | meaning |
|---|---|---|
| `REPRODUCED` | 155 | recomputed from the equations and the geometry, and agreeing to the printed resolution |
| `ARITHMETIC` | 93 | a deviation, a ratio, a dispersion or a count re-formed from the printed columns and found consistent with them |
| `NOT-FROM-MS` | 1 | an output of the coupled model; not recomputable from the manuscript |
| `SKIPPED (cost)` | 2 | the two longest tail sums of Table 5; the same constant is confirmed at the shorter truncations |

## How the two verdicts differ

`REPRODUCED` means the number was computed again, from the equations, and matched.

`ARITHMETIC` means the number is a function of other printed numbers, and the check is that it is
consistent with them. Those checks are made as **interval** tests, not as tolerance comparisons:
each printed input is taken to stand for the interval its own rounding leaves, the function is
evaluated over the corners of that box, and the printed result must lie within the printed
resolution of that interval. This is the correct test for a table whose deviations are formed at
full precision and rounded once at display (Section 5.3 of the manuscript): re-forming a deviation
from two rounded columns does not in general return the printed deviation, and the interval test is
what distinguishes that expected behaviour from an error.

## What was checked, block by block

| block | checks | what was recomputed |
|---|---|---|
| constants | 8 | (2μ₀L/π)ln 10, (2μ₀L/π)ln 2, 3μ₀L/(πd²), (3/4)·3μ₀L/(πd²), ln π + γ, U_a, U_m, the tail onset 1/U_a |
| Proposition 1, Table 8(a) | 11 | the four-argument closed form (29), and the direct summation at N = 10³…10⁷ |
| Table 8(b) | 14 | the diagonal by summation and by (30), and the increment per decade |
| Table 8(c) | 7 | the diagonal on the exact two-layer kernel, and the constant offset between the two kernels |
| Table 4 | 10 | the rank of Y and the null-space dimension at five truncations, from rank Y = min(2N, M_s − 1) |
| Table 5 | 8 | the residual tail of the hat series and its product by N², against 3μ₀L/(πd²) |
| Table 6(b) | 68 | both self-energy columns at twelve truncations, their increments, and the deviations from the two parameter-free constants |
| Table 6(a), 7(a) | 20 | the divergent bracket at every tiling and truncation, its growth over the released sweep, and its stationarity under the lock |
| Table 7(b) and continued | 22 | the twelve published dispersions, over six tilings and over the five finest, and the identity L_d = L_a − M |
| Table 7(d) | 8 | the eight deviations of the two bases against the reference |
| Table 12 | 52 | the fifty-one deviations, and the identity L_d = L_a − M on all four columns |
| Sections 7.4 and 8 | 2 | the eight–nine count over the seventeen quantities, and the thirteen-of-seventeen count |
| Section 7.1 | 9 | Carter's factor in both forms: g′, τ_s, x, γ, k_C, and the 0.16 % between them |
| Section 6.6 | 4 | the slotting frequency, the magnet conductivity, the skin depth and its ratio to the magnet thickness |
| Section 6.9 | 4 | the cost ratios: 0.84 %, 108 ms, ×115, 4.16 ms, and ×12.7 |

## The one value that is not reproducible from the manuscript alone

**Table 6(a), the two `B_g1` columns.** The fundamental of the air-gap flux density under a released
truncation is an output of the coupled model: it needs the reluctance network and the Newton solve,
neither of which is in this verification. What *is* checked on that panel is everything the
manuscript's own equations fix — the divergent bracket at each of the twelve truncations, the
increments, the ratios of successive increments, and the end-to-end drifts of 0.703 % and 0.023 %.

## Two findings

**One digit was wrong and has been corrected in the manuscript.** Table 12, row *MMF per magnet*,
column *mesh n_sh = 2*, printed a deviation of −1.32 %. The archived transcript
`outputs/MEC_BLDC/A1_table7_out.txt` carries the two full-precision values it is formed from,
759.32903 A against 769.44526 A, which give −1.31471 %. The printed figure does not follow from
them under any rounding, while the two neighbouring columns of the same row, −1.01 % and −0.85 %,
reproduce exactly. The manuscript now prints −1.31 %. Neither the eight–nine count of Section 7.4
nor the thirteen-of-seventeen count of Section 8 changes.

**The deviation column of Table 6(b) is not re-formable from the printed increment column.** At the
last doublings the increment is printed to seven significant figures, which is no longer enough to
re-form a deviation of order 10⁻⁵ %. The manuscript forms it at full precision, as Section 5.3
requires. The check here recomputes the self-energy at double precision and forms the increment and
the deviation from that, which reproduces the printed column. This is a property of the printed
table, not a defect: it is recorded so that a reader who tries the shorter route knows why it fails.

## What this does not establish

Agreement between two implementations of the same equations is not a proof that the equations
describe the machine, and it is not a validation of the coupled model. It establishes that the
closed forms of Section 3 are correctly stated, that the tables derived from them are correctly
computed, and that the arithmetic of the comparison tables is sound. The limits of the comparison
itself are in Section 7.4 of the manuscript and in `docs/OPEN_POINTS.md`.
