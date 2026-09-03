# Boundary condensation of a two-layer Dirichlet-to-Neumann operator

Reproduction archive for

> **Boundary Condensation of a Two-Layer Dirichlet-to-Neumann Operator: Trace Conformity and
> Air-Gap Mesh Reduction.**
> I. Laouar, A. Boukadoum, N. Mezhoud.
> Electrotechnical Laboratory of Skikda (LES), University 20 August 1955 Skikda, Algeria.
> Submitted to *Engineering Analysis with Boundary Elements*.

[![License](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Data](https://img.shields.io/badge/data-CC--BY--4.0-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-not%20yet%20minted-lightgrey.svg)](#citation)

This repository supports that one paper. It contains the MATLAB chain for the 750 W permanent-magnet
machine, the dated transcripts from which every published value is read, the finite-element exports
used as the reference, and a second, independent implementation of the paper's own equations.

---

## Start here

| I want to | Open |
|---|---|
| know what this archive can and cannot regenerate | [`REPRODUCIBILITY_STATUS.md`](REPRODUCIBILITY_STATUS.md) |
| check the closed forms without MATLAB | [`validation/INDEPENDENT_NUMERICAL_VERIFICATION.md`](validation/INDEPENDENT_NUMERICAL_VERIFICATION.md) |
| know where a published number comes from | [`docs/PROVENANCE.md`](docs/PROVENANCE.md) |
| re-run something | [`docs/REPRODUCING.md`](docs/REPRODUCING.md) |
| know what the paper does **not** establish | [`docs/OPEN_POINTS.md`](docs/OPEN_POINTS.md) |
| know which numbers moved, and why | [`CHANGELOG.md`](CHANGELOG.md) |

## Two tests, and neither needs a finite-element licence

**Without MATLAB — the equations.** A second implementation, in Python, of the closed forms of
Section 3 and of the tables that follow from them. It shares no code path with the MATLAB chain.

```bash
cd validation/independent
python3 verify_all.py
```

Expected: `checks: 251   agreeing: 248   mismatches: 0` and `VERDICT: PASS`, in about thirteen
seconds. Python 3.8 or later, standard library only.

**With MATLAB — the chain.** Reproduces a published quantity and checks itself.

```matlab
cd code
SET_REFERENCE_PATH          % rewrites the archived absolute paths, once, reversibly
cd MEC_BLDC
RUN_R5_NSH
```

Expected: `B_g1 = 1.07787 T` at `n_sh = 1` and `1.07908 T` at `n_sh = 2`, unknown counts 5400 and
6480, and the line `GARDE PASSEE`. Those two values are the first row of Table 12.

## Requirements

| | |
|---|---|
| Python | 3.8 or later, standard library only — for `validation/independent/` |
| MATLAB | R2024a; no toolbox beyond base MATLAB is required by the published chains |
| Reference solver | ANSYS Electronics Desktop 2023 R1, Maxwell 2-D, build 2023.1.0 — **not needed**: its exports are shipped under `reference/` as `.tab` tables |
| Disk | ≈ 7 MB |

**Runtime.** One end-to-end figure in this archive is reproducible, and it is the only one the paper
presents as such: the archived entry point runs in a median **18.91 s** over five repetitions on an
Intel Core i5-11400H at 2.70 GHz with 16 GB of memory (minimum 18.63 s, maximum 19.36 s, standard
deviation 0.29 s). Two further figures appear in the record, 143 s and 133 s, and **neither is
reproduced by any archived chain**; they are kept unadjusted and are not claimed as runtimes of this
chain. See `REPRODUCIBILITY_STATUS.md`.

The refined network of Section 6.9 is the polar mesh at `M_s = 540` with eleven interior layers plus
the ring of bore nodes: 540 × 12 = 6480 unknowns.

## Layout

| path | contents |
|---|---|
| `code/MEC_BLDC/` | the MATLAB chain for the 750 W permanent-magnet machine |
| `code/revision/` | the revision scripts, the manifest, and the two control checkers |
| `code/article/` | the figure-building scripts and a **previous** version of the manuscript source |
| `code/SET_REFERENCE_PATH.m` | rewrites the archived absolute paths to this checkout |
| `validation/independent/` | the second implementation, in Python, and its results |
| `scripts/` | the entry points, and what each one runs |
| `environment/` | what has to be installed, and what does not |
| `outputs/MEC_BLDC/` | execution transcripts (`diary`), at full precision |
| `outputs/revision/` | transcripts of the revision blocks B1 to B12 |
| `outputs/annulus_reference/` | transcripts on a thin annulus, kept for the reason below |
| `reference/ANSYS_750W/` | the finite-element exports used as the reference |
| `manuscript/` | the submitted text |
| `docs/` | provenance, reproduction guide, open points, independent audit |
| `notes/` | delivery notes for each verification block (in French; working documents) |

## Why a thin-annulus transcript is in this archive

`outputs/annulus_reference/` holds transcripts produced on a thin annulus, with the dimensions of an
18.5 kW induction machine. They are kept because the closed-form divergence proved in Section 3 is
observable in the outputs only in the thin-annulus regime. **No value in the submitted paper is read
from them**, and the manifest records, for the one block that resembles a published table, why it
does not correspond to it. On the permanent-magnet machine the truncation is locked to the tiling by
the rank condition, so the bracket is stationary and the divergence cannot be exhibited on that chain
without releasing the truncation — which the paper does, and reports.

## How a number is traced

Every value in the paper comes from a transcript under `outputs/`, or is declared as not coming from
one. The project rule is **one quantity, one chain**: where two scripts produced the same quantity,
one was archived and the other withdrawn rather than kept alongside.

Each transcript opens with the configuration that produced it: machine, tiling `M_s`, radial layers
`n_sh`, surface basis, solver, fringing constant, and the path of the finite-element reference.

Reading rule: the authoritative block is the **last complete** block of a transcript — complete
meaning its internal consistency has been checked, not merely that it comes last.

Deviations in the paper are formed on full-precision values and rounded once, at display. Re-forming
a deviation from two rounded printed columns will not in general reproduce the printed deviation, and
that is expected; `validation/independent/verify_all.py` tests exactly that, as an interval check.

## Verification blocks

Each block in `notes/` carries a guard: a check that would contradict its own result if the result
were wrong. A block whose guard was not executed was not accepted.

- `R4_NOTE.md` — unlocks the truncation from the tiling; the drift returns at the predicted slope.
- `R5_NOTE.md` — settles which shoe layering produced the published column.
- `R6_NOTE.md` — examines the saturation map and refutes a proposed explanation of the sweep deficit.
- `R7_NOTE.md` — recounts the twenty-nine compared quantities, and reports a defect found in doing so.
- `R9_NOTE.md` — archive assembly.
- `G_NOTE.md` — closure of 26 August 2026: six guards, none failed.

`docs/AUDIT_INDEPENDANT.md` reports an independent audit of the archive; `validation/` reports the
independent recomputation of the equations.

## Citation

**No DOI has been minted, and none has been invented.** The paper's data-availability statement says
the persistent identifier is registered on acceptance and printed in the published version. Once the
deposit is made, write the DOI — unchanged — in four places: this file, `CITATION.cff`,
`.zenodo.json`, and the manuscript.

```bibtex
@software{laouar_dtn_pmsm_archive,
  author  = {Laouar, Idris and Boukadoum, Ahcene and Mezhoud, Nabil},
  title   = {Reproduction archive for the boundary condensation of a two-layer
             Dirichlet-to-Neumann operator},
  year    = {2026},
  version = {1.2.0},
  note    = {DOI to be registered on acceptance},
  url     = {https://github.com/il1996/mec-dtn-airgap-PMSM}
}
```

## Licence

Code under MIT, data and documents under CC BY 4.0. See [`LICENSE`](LICENSE).
**The `ACTION REQUIRED` notice in that file must be resolved before the repository is made public.**
