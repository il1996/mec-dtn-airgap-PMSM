# A Well-Posed Dirichlet-to-Neumann Air-Gap Operator for Magnetic Equivalent Circuits

Reproduction repository for the paper

> **A Well-Posed Dirichlet-to-Neumann Air-Gap Operator for Magnetic Equivalent
> Circuits: Trace Conformity, and Validation on a 15-Slot/14-Pole PMSM.**
> I. Laouar, A. Boukadoum, N. Mezhoud.
> Electrotechnical Laboratory of Skikda (LES), University 20 August 1955 Skikda, Algeria.

[![License](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Data](https://img.shields.io/badge/data-CC--BY--4.0-blue.svg)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-not%20yet%20minted-lightgrey.svg)](#citation)

This repository supports **one paper only**. The companion paper on the doubly
slotted annulus operator is deposited separately; neither archive depends on
the other to be read or re-run.

---

## Start here

Two files answer most questions a reader will have:

| I want to | Open |
|---|---|
| know where a published number comes from | [`docs/PROVENANCE.md`](docs/PROVENANCE.md) |
| re-run something | [`docs/REPRODUCING.md`](docs/REPRODUCING.md) |
| know what the paper does **not** establish | [`docs/OPEN_POINTS.md`](docs/OPEN_POINTS.md) |
| know which numbers moved, and why | [`CHANGELOG.md`](CHANGELOG.md) |

## The two-minute test

Run this before anything else. It reproduces a published quantity and checks
itself.

```matlab
cd code
SET_REFERENCE_PATH          % rewrites the archived absolute paths, once, reversibly
cd MEC_BLDC
RUN_R5_NSH
```

Expected: `B_g1 = 1.07787 T` at `n_sh = 1` and `1.07908 T` at `n_sh = 2`;
unknown counts 5400 and 6480. The script ends with `GARDE PASSEE`.

If that passes, the archive works on your machine.

## Requirements

| | |
|---|---|
| MATLAB | R2024a (no toolbox beyond base MATLAB is required by the published chains) |
| Reference solver | ANSYS Electronics Desktop 2023 R1, Maxwell 2-D (build 2023.1.0) — **not needed to re-run the network**; its exports are included under `reference/` |
| Disk | ≈ 5 MB |
| Runtime | the refined network is 10 800 unknowns and runs in ≈ 133 s on an Intel Core i5-11400H at 2.70 GHz with 16 GB of memory |

Nothing here needs a licence for the finite-element solver. The reference is
shipped as exported `.tab` tables, so every comparison in the paper can be
re-formed without re-solving it.

## Layout

| path | contents |
|---|---|
| `code/MEC_BLDC/` | the MATLAB chain for the 750 W permanent-magnet machine |
| `code/article/` | manuscript source, bibliography, and the scripts that build the figures |
| `code/SET_REFERENCE_PATH.m` | rewrites the archived absolute paths to this checkout |
| `outputs/MEC_BLDC/` | execution transcripts (`diary`), at full precision |
| `outputs/annulus_reference/` | transcripts on a thin induction-machine annulus — see below |
| `reference/ANSYS_750W/` | the finite-element exports used as the reference |
| `docs/` | provenance, reproduction guide, open points, independent audit |
| `notes/` | delivery notes for each verification block (in French; working documents) |
| `CHANGELOG.md` | what changed at each release, and which numbers moved |

## Why an induction-machine annulus is in a permanent-magnet paper

`outputs/annulus_reference/` holds transcripts produced on a thin annulus taken
with the dimensions of an 18.5 kW induction machine. This is deliberate and is
stated in the paper: the closed-form divergence proved there is observable only
in the thin-annulus regime. On the permanent-magnet machine the truncation is
locked to the tiling by the rank condition, so the bracket is stationary and the
divergence cannot be exhibited on that chain without releasing the truncation —
which the paper also does, and reports.

## How a number is traced

Every value in the paper comes from a transcript under `outputs/`. The project
rule is **one quantity, one chain**: where two scripts produced the same
quantity, one was archived and the other withdrawn rather than kept alongside.

Each transcript opens with the configuration that produced it: machine, tiling
`M_s`, radial layers `n_sh`, surface basis, solver, fringing constant, and the
path of the finite-element reference used.

Reading rule: the authoritative block is the **last complete** block of a
transcript — complete meaning its internal consistency has been checked, not
merely that it comes last. Position alone is not a criterion.

## Verification blocks

Each block in `notes/` carries a guard: a check that would contradict its own
result if the result were wrong. A block whose guard was not executed was not
accepted.

- `R4_NOTE.md` — unlocks the truncation from the tiling; the drift returns, at
  the predicted slope to machine precision.
- `R5_NOTE.md` — settles which shoe layering produced the published column.
- `R6_NOTE.md` — examines the saturation map and refutes a proposed explanation
  of the speed-sweep deficit.
- `R7_NOTE.md` — recounts the twenty-nine compared quantities, and reports a
  defect found while recounting.
- `R9_NOTE.md` — archive assembly.
- `G_NOTE.md` — **closure, 26 August 2026**: settles the convention on the
  Euler constant in the divergent bracket, reads the four bracket widths back at
  full precision, publishes the complete matrix of the reluctance factor,
  converts the two efficiency rows to points and recounts everything that
  depends on them. Six guards, none failed.

`docs/AUDIT_INDEPENDANT.md` reports an independent audit of the archive.

## Citation

**No DOI has been minted, and none has been invented.** The manuscript carries
no unresolvable identifier either: its data-availability statement says the
persistent identifier is registered on acceptance and printed in the published
version.

Once the deposit is made, write the DOI — unchanged — in four places: this
file, `CITATION.cff`, `.zenodo.json`, and the manuscript. Until then, cite the
repository by its URL:

```bibtex
@software{laouar_dtn_pmsm_archive,
  author  = {Laouar, Idris and Boukadoum, Ahcene and Mezhoud, Nabil},
  title   = {Reproduction archive for a well-posed Dirichlet-to-Neumann
             air-gap operator for magnetic equivalent circuits},
  year    = {2026},
  note    = {DOI to be registered on acceptance},
  url     = {https://github.com/USER/REPO}
}
```

## Licence

Code under MIT, data and documents under CC BY 4.0. See [`LICENSE`](LICENSE).
**The copyright line must be completed before the repository is made public.**
