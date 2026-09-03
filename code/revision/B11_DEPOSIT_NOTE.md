# B11 — Deposit: what was audited, what was prepared, and what only you can do

**Date:** 2026-09-03  ·  **Task:** B11 of the Deliverable 3 specification.

This note is the deliverable for B11. Three of its four items are audits and
proposed replacements, which are below. One item — obtaining the DOI — cannot be
done from here, and the reason is stated rather than worked around.

---

## 1. Zenodo deposit and the DOI — **blocked, and it is yours to unblock**

The specification asks for a Zenodo deposit with the **DOI obtained at
submission** rather than "registered on acceptance".

I have not created the deposit. Minting a DOI means creating or using an account
on an external service and publishing this archive under your name; that is an
outward-facing, effectively irreversible act on your behalf, and it needs your
hand on it, not mine. The metadata it needs is prepared below so the deposit is
a short operation when you do it.

**What to do, in order.**

1. Reserve the DOI on Zenodo *before* uploading (Zenodo's "Reserve DOI" button on
   a draft deposit). This is what makes the identifier citable at submission.
2. Upload the archive.
3. Write the DOI, unchanged, in **four** places — `README.md`, `CITATION.cff`,
   `.zenodo.json`, and the manuscript's data-availability statement.

---

## 2. `.zenodo.json` — **wrong paper**

The archived metadata still describes the *previous* manuscript:

```
"title": "Reproduction archive for a well-posed Dirichlet-to-Neumann air-gap
          operator for magnetic equivalent circuits (Paper I)"
"description": "... for the paper <em>A Well-Posed Dirichlet-to-Neumann Air-Gap
          Operator for Magnetic Equivalent Circuits: Trace Conformity, and
          Validation on a 15-Slot/14-Pole PMSM</em> ..."
```

The submitted manuscript is **“Boundary Condensation of a Two-Layer
Dirichlet-to-Neumann Operator: Trace Conformity and Air-Gap Mesh Reduction.”**
An editor comparing the deposit to the submission would see two different papers.

`CITATION.cff` carries the same old title, already at `version: "1.1.0"`, with the
DOI line correctly commented out rather than invented.

**I have not edited either file.** They are published repository metadata, the
DOI cannot be filled yet, and changing the title is a decision about how you
present the deposit. Both need one edit each, at the same time as the DOI.

---

## 3. Data availability — the stray full stop, and more

Current text in the manuscript:

> "...deposited in a public archive whose persistent identifier is registered on
> acceptance.**;** pending that, at https://github.com/il1996/mec-dtn-airgap-PMSM
> (v1.0.0)."

Three defects, not one:

| # | Defect | Note |
|---|---|---|
| 1 | `acceptance.;` — a full stop inside the sentence | the one the specification names |
| 2 | the URL carries a personal GitHub username | `il1996` |
| 3 | `(v1.0.0)` contradicts `CITATION.cff`, which says `version: "1.1.0"` | a version mismatch between the manuscript and the archive |

Proposed replacement, to be applied **in the editorial pass, not by me** (this
task set is computational; the manuscript is not edited here):

> The programs generating every figure and table, the diagnostic scripts for the
> hypotheses rejected in Sections 7.2 and 7.3, the manifest binding each
> published quantity to its chain, the uncertainty budget of Section 7.2 and the
> full table of the twenty-nine compared quantities are deposited at
> `https://doi.org/10.5281/zenodo.XXXXXXX` (version 1.1.0).

---

## 4. Redistributability — **audited, and it is clean, with one caveat**

**What is in `reference/`:** 38 files, *all* `.tab`. They are plain tab-separated
text exports — the first line of one reads `"$OSD [mm]"  "k_l []"` followed by
numbers. Nothing proprietary to ANSYS is redistributed: there is **no** `.aedt`,
`.aedtresults`, `.aedz`, `.asol`, `.lock`, `.pjt`, `.exe` or `.dll` anywhere in
the repository. The finite-element solver is not needed to re-run anything.

**What is deliberately *not* redistributed**, and should be said so explicitly:
five Word files present in the original project folder and absent from the
archive —

```
magnetostique(Armature-Field)/résultat.docx
magnetostique(Magnetic_loading)/résultat.docx
transitoire (a vide)/résultat.docx
transitoire (Back_emf)/résultats.docx
transitoire (en charge)/résultat.docx
```

These are working notes, not data any published value depends on. No script
reads them. Their absence is correct; it should simply be stated rather than
silent.

**The caveat, and it matters to a stranger.** `code/MEC_BLDC/machine_bldc.m`
line 108 hard-codes an absolute path *outside* the repository:

```matlab
M.FEA.dir = 'C:\Users\hp\Desktop\ANSYS résultat 750W';
```

This is deliberate and documented — `code/SET_REFERENCE_PATH.m` exists precisely
to rewrite it, reversibly, with a `.bak` beside each file — and the reasoning
(archive the scripts byte-identical to those that produced the numbers) is sound.
But **the archive does not work out of the box**: any script that reads the
finite-element reference fails until `SET_REFERENCE_PATH` has been run. The
README says to run it; the deposit description does not. It should.

---

## 5. Tag `v1.1.0` — **not created, and here is why**

The specification asks for a tag `v1.1.0` matching the revised manuscript.

`git -C mec-dtn-airgap-pmsm log` reports:

```
fatal: your current branch 'main' does not have any commits yet
```

Every file is **staged but never committed**. Creating the tag therefore means
first creating the repository's initial commit — which decides the authorship,
the message and the history of your archive. That is a call for you, not a side
effect of a verification task, so I have staged nothing further and committed
nothing.

When you want it:

```bash
git -C mec-dtn-airgap-pmsm add code/revision outputs/revision
git -C mec-dtn-airgap-pmsm commit -m "Reproduction archive v1.1.0"
git -C mec-dtn-airgap-pmsm tag -a v1.1.0 -m "Archive matching the revised manuscript"
```

Note that `v1.1.0` should be cut **after** the B1–B12 artefacts are in, since
those are part of what the revised manuscript will cite.

---

## Summary of B11

| Item | State |
|---|---|
| Zenodo deposit and DOI at submission | **blocked — needs your account**; metadata prepared |
| Tag `v1.1.0` | **not created** — the repository has no commits yet; command given |
| `Data availability` rewritten | proposed text above; three defects found, not one |
| Nothing non-redistributable is depended on | **verified** — 38 plain-text `.tab` exports, no proprietary artefacts; five working `.docx` correctly excluded; the out-of-repo absolute path is documented but should be flagged in the deposit description |

---

## Postscript — 3 September 2026

This note is kept as the dated record it was. Two of its statements have since been overtaken and
are corrected here rather than edited above:

- **“the repository has no commits yet”** was true when the note was written. The archive was
  committed and tagged on 3 September 2026; see `CHANGELOG.md`.
- **The tag was cut as `v1.2.0`, not `v1.1.0`.** The 1.1.0 closure of 26 August 2026 was never
  released on its own; the pre-submission audit of 3 September 2026 followed it, and the two were
  published together. `CHANGELOG.md` carries both entries.

The DOI is still not minted and still not invented. `docs/OPEN_POINTS.md` §6 remains open.
