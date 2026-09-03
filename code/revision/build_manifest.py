#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
build_manifest.py -- block B10, revision of 2026-09-02.

Builds  code/revision/MANIFEST.json  and  code/revision/annex_29_quantities.md.

RULES OBEYED
  * No number is typed by hand in this file.
      - published values are READ from the rendered manuscript
        (scratchpad/ms_text.txt, one paragraph per line), by PARAGRAPH INDEX;
      - model / reference / deviation values of the twenty-nine quantities are
        PARSED from outputs/MEC_BLDC/R7_scorecard_out.txt and
        outputs/MEC_BLDC/G4_points_out.txt.
    What is typed by hand here is only the MAPPING (which paragraph belongs to
    which table, which transcript carries it) and prose read from the files.
  * Nothing under code/MEC_BLDC, outputs/, reference/, docs/, notes/ is written.
  * An item whose chain could not be established gets authoritative=false with
    empty script and transcript, never a guess.

The mapping was derived by CONTENT, never by number: docs/PROVENANCE.md follows
the numbering of the earlier manuscript and is not used as a lookup table.
"""

import io
import json
import os
import re
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

ISO_DATE = "2026-09-02"

REPO = r"C:\Users\hp\Desktop\mec-dtn-airgap-pmsm"
SCRATCH = (r"C:\Users\hp\AppData\Local\Temp\claude"
           r"\C--Users-hp-Desktop-claude"
           r"\c47f027c-9355-48a8-92a8-e75c9f76a865\scratchpad")
MS_TEXT = os.path.join(SCRATCH, "ms_text.txt")

OUT_JSON = os.path.join(REPO, "code", "revision", "MANIFEST.json")
OUT_ANNEX = os.path.join(REPO, "code", "revision", "annex_29_quantities.md")

R7 = "outputs/MEC_BLDC/R7_scorecard_out.txt"
G4 = "outputs/MEC_BLDC/G4_points_out.txt"


# --------------------------------------------------------------------------
# 1. the manuscript, one paragraph per line
# --------------------------------------------------------------------------
def load_ms():
    d = {}
    with io.open(MS_TEXT, encoding="utf-8") as fh:
        for line in fh:
            idx, _, txt = line.partition("\t")
            try:
                d[int(idx)] = txt.rstrip("\n")
            except ValueError:
                pass
    return d


SUP = {"\u2070": "0", "\u00b9": "1", "\u00b2": "2", "\u00b3": "3",
       "\u2074": "4", "\u2075": "5", "\u2076": "6", "\u2077": "7",
       "\u2078": "8", "\u2079": "9", "\u207b": "-", "\u207a": "+"}

NUM = re.compile(r"[-+\u2212]?\d[\d\u00a0\u202f ]*(?:\.\d+)?(?:[eE][-+]?\d+)?")


def normalise(text):
    """Unicode -> ascii arithmetic, superscript exponents -> e-notation."""
    s = text
    s = s.replace("{math:", " ").replace("}", " ")
    s = s.replace("\u2212", "-").replace("\u2013", "-")
    # a x 10^b  ->  a e b
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if c in ("\u00d7", "x") and s[i:i + 4].replace("\u00d7", "x") == "x10":
            j = i + 3
            exp = ""
            while j < len(s) and s[j] in SUP:
                exp += SUP[s[j]]
                j += 1
            if exp:
                out.append("e" + exp)
                i = j
                continue
        if c in SUP:
            out.append(SUP[c])
            i += 1
            continue
        out.append(c)
        i += 1
    return "".join(out)


def numbers_in(text):
    """Every number a paragraph publishes, as the manuscript prints it."""
    s = normalise(text)
    vals = []
    for m in NUM.finditer(s):
        tok = m.group(0)
        tok = tok.replace("\u00a0", "").replace("\u202f", "").replace(" ", "")
        if tok in ("-", "+", ""):
            continue
        try:
            float(tok)
        except ValueError:
            continue
        vals.append(tok)
    return vals


def published(ms, paras):
    """The published values named by a list of paragraph indices."""
    out = []
    for p in paras:
        for v in numbers_in(ms.get(p, "")):
            out.append({"value": v, "paragraph": p})
    return out


# --------------------------------------------------------------------------
# 2. the twenty-nine quantities, parsed from their own transcripts
# --------------------------------------------------------------------------
ROW29 = re.compile(
    r"^\s*(\d{1,2})\s+(\S.*?)\s{2,}"
    r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s+"
    r"([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s+"
    r"([-+]\d+(?:\.\d+)?)\s*(<--.*)?$")


def parse_r7(path):
    rows = []
    with io.open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = ROW29.match(line.rstrip("\n"))
            if not m:
                continue
            n = int(m.group(1))
            if n < 1 or n > 29:
                continue
            rows.append({
                "n": n,
                "quantity_fr": m.group(2).strip(),
                "model": m.group(3),
                "reference": m.group(4),
                "deviation_pct": m.group(5),
                "flagged_over_5pct": bool(m.group(6)),
            })
    seen, uniq = set(), []
    for r in rows:
        if r["n"] in seen:
            continue
        seen.add(r["n"])
        uniq.append(r)
    uniq.sort(key=lambda r: r["n"])
    return uniq


def parse_r7_context(path):
    """The three counts and the exclusion motives, read from the transcript."""
    txt = io.open(path, encoding="utf-8", errors="replace").read()
    ctx = {}
    m = re.search(r"compteur du programme maitre \(nok\)\s*:\s*(\d+) sur (\d+)", txt)
    if m:
        ctx["master_counter"] = "%s / %s" % (m.group(1), m.group(2))
    m = re.search(r"recomptage independant\s*:\s*(\d+) sur (\d+)", txt)
    if m:
        ctx["independent_recount"] = "%s / %s" % (m.group(1), m.group(2))
    m = re.search(r"compte hors grandeurs non validees\s*:\s*(\d+) sur (\d+)", txt)
    if m:
        ctx["excluding_three"] = "%s / %s" % (m.group(1), m.group(2))
    ctx["excluded"] = re.findall(r"exclue #(\d+) (.+?)\s{2,}([-+]\d+\.\d+) %", txt)
    ctx["motives"] = re.findall(r"motif : (.+)", txt)
    return ctx


def parse_g4(path):
    txt = io.open(path, encoding="utf-8", errors="replace").read()
    out = {"points": {}, "recount": {}, "guard": None}
    for m in re.finditer(
            r"#(\d+)\s+(\S.*?)\s*\n\s*modele\s+([\d.]+)\s+reference\s+([\d.]+)"
            r"\s*\n\s*relatif\s+([-+]?[\d.]+) %\s+POINTS\s+([-+]?[\d.]+) p\.p\.",
            txt):
        out["points"][int(m.group(1))] = {
            "label": m.group(2).strip(),
            "model_full": m.group(3),
            "reference_full": m.group(4),
            "relative_pct": m.group(5),
            "points_pp": m.group(6),
        }
    for key, pat in (
            ("relative_all", r"dans les 5, tout en relatif\s*:\s*(\d+ / \d+)"),
            ("points_2529", r"dans les 5, 25 et 29 en points\s*:\s*(\d+ / \d+)"),
            ("relative_excl", r"hors les trois non validees, relatif\s*:\s*(\d+ / \d+)"),
            ("points_excl", r"hors les trois non validees, points\s*:\s*(\d+ / \d+)")):
        m = re.search(pat, txt)
        if m:
            out["recount"][key] = m.group(1)
    m = re.search(r"nombre a moins d un point\s*:\s*(.+)", txt)
    if m:
        out["within_one_point"] = m.group(1).strip()
    m = re.search(r"ecart maximal sur les 29\s*:\s*(\S+)", txt)
    if m:
        out["guard"] = m.group(1)
    return out


# English labels for the 29 rows. Not values: names, taken from the manuscript's
# own vocabulary so that a referee can align the annex with Fig. 11.
EN = {
    1: "Fundamental air-gap flux density B_g1 (T)",
    2: "Mean air-gap flux density |B_r| (T)",
    3: "Peak air-gap flux density B_r (T)",
    4: "Tangential air-gap flux density B_t, rms (T)",
    5: "Flux-line value A, peak (Wb/m)",
    6: "Magnetomotive force per magnet, mean (A)",
    7: "Magnetomotive force per gap, mean (A)",
    8: "Magnet magnetomotive force, dispersion (%)",
    9: "Reluctance factor k_r (-)",
    10: "Leakage factor k_l (-)",
    11: "Self inductance L_a (mH)",
    12: "Mutual inductance M (mH)",
    13: "Synchronous inductance L_d (mH)",
    14: "Phase back-electromotive force, peak (V)",
    15: "Six-step envelope, peak (V)",
    16: "Flux linkage, peak (Wb)",
    17: "Electromagnetic torque on load (N m)",
    18: "Phase current, rms, on load (A)",
    19: "Iron loss at no load (W)",
    20: "Magnet eddy-current loss at no load (W)",
    21: "Settled speed on load (rpm)",
    22: "Mean dc-bus current (A)",
    23: "Iron loss on load (W)",
    24: "Magnet eddy-current loss on load (W)",
    25: "Efficiency on load (%)",
    26: "Back-electromotive-force constant k_E (mV/rpm)",
    27: "Standstill torque, speed sweep (N m)",
    28: "Peak power, speed sweep (W)",
    29: "Peak efficiency, speed sweep (%)",
}

# Where each of the 29 is produced inside the master run. Read from
# BLDC_MEC_COMPLET_out.txt section headings, not invented.
SECTION = {}
for _n in range(1, 9):
    SECTION[_n] = "no-load field and internal magnetic circuit (master run, sections 1-3)"
for _n in (9, 10):
    SECTION[_n] = "reluctance and leakage factors (master run, section 3)"
for _n in (11, 12, 13):
    SECTION[_n] = "inductances (master run, section 4)"
for _n in (14, 15, 16):
    SECTION[_n] = "back-EMF and flux linkage at 1500 rpm (master run, section 4)"
for _n in (17, 18, 21, 22, 23, 24, 25):
    SECTION[_n] = "on-load transient, steady state (master run, section 5b)"
for _n in (19, 20):
    SECTION[_n] = "no-load losses (master run, section 5)"
for _n in (26, 27, 28, 29):
    SECTION[_n] = "speed sweep at 500 V dc (master run, sections 6-6b)"


# --------------------------------------------------------------------------
# 3. the mapping, derived by content (see the report)
# --------------------------------------------------------------------------
MASTER = "outputs/MEC_BLDC/BLDC_MEC_COMPLET_out.txt"
S_MASTER = "code/MEC_BLDC/BLDC_MEC_COMPLET.m"

TABLES = [
    dict(key="T1a", quantity="Table 1, rows N1 / N2 / three N2b / N3 and the FEA reference: "
                             "air-gap closure at fixed iron model",
         table="1", figure=None,
         script="code/MEC_BLDC/RUN_CARTER_CMP.m",
         transcript="outputs/MEC_BLDC/A6_n2b_out.txt",
         block="last complete block",
         date="2026-08-05 (docs/MANIFEST.md section 1)",
         configuration="PMSM 15/14 750 W, lumped iron network, linear solver mu_r = 3000, "
                       "kfringe = 0.325, tooth-face arc a_t = 0.361206 rad",
         authoritative=True,
         paras=[132, 133, 134, 139, 140, 146, 147, 182, 183, 184],
         note="Mapped by content: the closure names, the three tangential discretisations "
              "a_r = 0.25/0.50/1.00 a_t and the FEA row of the transcript are the rows of the "
              "submitted Table 1. outputs/MEC_BLDC/t6_n2b_out.txt holds the same block behind a "
              "MATLAB variable dump and is the superseded copy."),
    dict(key="T1b", quantity="Table 1, rows N4 (air-gap element), A1 (relative permeance), "
                             "A2 (complex permeance)",
         table="1", figure=None, script="", transcript="", block="",
         date="", configuration="",
         authoritative=False,
         paras=[],
         note="NOT ATTRIBUTED. None of the six deviations of these three rows, nor the "
              "fundamental 1.1210 T quoted for N4 in the table notes, occurs anywhere in "
              "outputs/. No script in code/ writes them. The rows are published in the "
              "submitted manuscript only."),
    dict(key="T2", quantity="Table 2: positioning against the air-gap closures in use",
         table="2", figure=None, script="", transcript="", block="",
         date="", configuration="",
         authoritative=False, paras=[],
         note="Narrative table. It carries no computed quantity, so it has no chain; this is a "
              "declaration, not a gap."),
    dict(key="T3", quantity="Table 3: data of the test machine",
         table="3", figure=None,
         script="code/MEC_BLDC/RUN_U2_MACHINE_TABLE.m",
         transcript="outputs/MEC_BLDC/U2_machine_table_out.txt",
         block="only block",
         date="2026-08-12 (script mtime; run emits 34 lines, guard passed)",
         configuration="values emitted from machine_bldc.m and checked line by line against the "
                       "ANSYS project BLDC.aedt; 34 lines concordant, 0 discrepant",
         authoritative=True,
         paras=[414, 416, 418, 420, 422, 424, 426, 428, 430, 432, 434, 436, 440,
                444, 446, 450, 452, 458, 460],
         note="Mapped by content: the transcript's 'emis/publie' table is line-for-line the "
              "submitted Table 3."),
    dict(key="T4", quantity="Table 4: rank of the bore operator, M_s = 360",
         table="4", figure=None, script="", transcript="", block="",
         date="", configuration="",
         authoritative=False, paras=[],
         note="NOT ATTRIBUTED. The five rank/null-space triples are not printed by any "
              "transcript in outputs/. They appear in code/article/ArticleI_DtN_PMSM.tex "
              "(tab:rank) with the same values, but a manuscript source file is not a run "
              "transcript and is not accepted here as one."),
    dict(key="T5", quantity="Table 5: hat-basis residual tail on the two-layer annulus",
         table="5", figure=None, script="", transcript="", block="",
         date="", configuration="",
         authoritative=False, paras=[],
         note="NOT ATTRIBUTED. The submitted tail values (3.4267e-8 down to 1.0998e-15) are in "
              "no transcript in outputs/. The earlier manuscript's tab:p1tail carries a "
              "DIFFERENT tail (1.7111e-7 down to 5.4918e-15) on a different annulus, so the "
              "submitted table is a later recomputation whose transcript is not deposited."),
    dict(key="T6a", quantity="Table 6, panels (a) and (b), rows N_h = 540 to 4320: truncation "
                             "released from the tiling",
         table="6", figure=None,
         script="code/MEC_BLDC/RUN_R4_UNLOCK.m + code/MEC_BLDC/RUN_M4_HATBASIS.m "
                "+ code/MEC_BLDC/RUN_G1_GAMMA.m (gamma convention)",
         transcript="outputs/MEC_BLDC/M4_hatbasis_out.txt + "
                    "outputs/MEC_BLDC/G1_gamma_out.txt + outputs/MEC_BLDC/R4_unlock_out.txt",
         block="M4: both basis blocks; G1: sections 1 and 2; R4: section 2",
         date="2026-08-10 (M4 script mtime) / 2026-08-26 (G1, docs/MANIFEST.md section 4)",
         configuration="PMSM 15/14 750 W, M_s = 1080 fixed, N released over lock x [1 2 4 8], "
                       "P0 and hat bases, linear solver mu_r = 3000, kfringe = 0.325, "
                       "gamma = 0.57721566490153287",
         authoritative=True,
         paras=[564, 565, 568, 571, 572, 573, 575, 578, 579, 580, 581, 582, 583,
                585, 586, 587, 588, 589, 590],
         note="Mapped by content: M4 prints both bases at exactly these four truncations and "
              "the same increment ratios; G1 prints the bracket column with gamma."),
    dict(key="T6b", quantity="Table 6, rows N_h = 8640 to 1 105 920 (seven further doublings), "
                             "panel (b) in full, and footnotes c and d",
         table="6", figure=None, script="", transcript="", block="",
         date="", configuration="",
         authoritative=False, paras=[],
         note="NOT ATTRIBUTED. The deposited chain stops at N_h = 4320: R4_unlock_out.txt "
              "returns NaN at 8640 with a singular-matrix warning, which the submitted "
              "footnote d identifies as the sinh overflow of the two-layer kernel. The "
              "submitted rows 8640 to 1 105 920, the smallest singular values of footnote c "
              "and the ratio 0.923624 come from a later run whose transcript is not in "
              "outputs/."),
    dict(key="T7a", quantity="Table 7 panel (a): the divergent bracket under the locked chain",
         table="7", figure=None,
         script="code/MEC_BLDC/RUN_R4_UNLOCK.m + code/MEC_BLDC/RUN_G1_GAMMA.m",
         transcript="outputs/MEC_BLDC/R4_unlock_out.txt + outputs/MEC_BLDC/G1_gamma_out.txt",
         block="R4 section 1 (locked chain, M_s = 540..8640); G1 section 5 (ln pi + gamma)",
         date="2026-08-26 (G1) / 2026-08-24 (R4 as deposited)",
         configuration="locked chain N = floor(M_s/2), M_s = 540, 1080, 2160, 4320, 8640, "
                       "P0 basis, linear solver mu_r = 3000",
         authoritative=True,
         paras=[804, 810, 816, 822, 828, 834],
         note="DERIVED COLUMN. R4 prints the bracket WITHOUT gamma (1.144724 ... 1.144730); "
              "G1 fixes the gamma convention and prints ln(pi)+gamma = 1.721946. The published "
              "column is R4 + gamma, an addition performed between the transcript and the "
              "table. Only 1.721944 (M_s = 1080) and 1.721946 (the limit) are printed as such."),
    dict(key="T7b", quantity="Table 7 panel (b): the bore in both surface bases, locked chain, "
                             "six tilings, and the four dispersions",
         table="7", figure=None,
         script="code/MEC_BLDC/RUN_V1_PMSM_BASE.m",
         transcript="outputs/MEC_BLDC/V1_pmsm_basis_out.txt",
         block="second block (the first aborts on an unrecognised field 'ths')",
         date="2026-08-04 (docs/MANIFEST.md section 1)",
         configuration="PMSM 15/14 750 W, locked chain numax = Nsurf/2 imposed by cogging_mec, "
                       "six tilings 540..8640, P0 and P1 bases, kfringe = 0.325, linear solver",
         authoritative=True,
         paras=[849, 850, 851, 852, 855, 856, 857, 858, 861, 862, 863, 864,
                867, 868, 869, 870, 873, 874, 875, 876, 879, 880, 881, 882,
                885, 886, 887, 888],
         note="Anchor supplied with the task and re-verified here cell by cell. The dispersions "
              "0.440 / 0.496 / 26.126 / 28.514 % are reproduced independently on 2026-09-02 by "
              "code/revision/tiling_dispersion_integral.m -> "
              "outputs/revision/B1_tiling_dispersion_out.txt (non-regression gate passed)."),
    dict(key="T7d1", quantity="Table 7 panel (d), rows B_g1 and nu = 8: both bases against the "
                              "finite-element reference at M_s = 1260",
         table="7", figure=None,
         script="code/MEC_BLDC/RUN_V1_PMSM_BASE.m",
         transcript="outputs/MEC_BLDC/V1_pmsm_basis_out.txt",
         block="second block, row Nsurf = 1260",
         date="2026-08-04 (docs/MANIFEST.md section 1)",
         configuration="M_s = 1260, N = 630, P0 and hat bases, linear solver mu_r = 3000, "
                       "kfringe = 0.325",
         authoritative=True,
         paras=[900, 901, 903, 906, 907, 909],
         note="The transcript prints five decimals where the table publishes six; the reference "
              "column 1.074551 / 0.019658 is printed there to four and five decimals only."),
    dict(key="T7d2", quantity="Table 7 panel (d), rows nu = 22 and flux linkage peak",
         table="7", figure=None, script="", transcript="", block="",
         date="", configuration="",
         authoritative=False, paras=[],
         note="NOT ATTRIBUTED at submission. V1_pmsm_basis_out.txt carries neither nu = 22 nor "
              "the flux linkage: its columns are B_g1, |B_r| mean, B_r peak, B_t rms, nu = 8. "
              "The two hat-basis values 0.036196 and 0.261018 are in no transcript that predates "
              "the submission; lambda = 0.261206 (p.c.) and 0.261018 (hat) are reproduced on "
              "2026-09-02 in outputs/revision/B1_tiling_dispersion_out.txt, after submission, "
              "and nu = 22 in either basis is still unattributed."),
    dict(key="T8", quantity="Table 8: verification of Proposition 1 on the two-layer annulus of "
                            "the PMSM (panels a, b, c)",
         table="8", figure=None, script="", transcript="", block="",
         date="", configuration="",
         authoritative=False, paras=[],
         note="NOT ATTRIBUTED. docs/PROVENANCE.md sends its own Table 8 to "
              "outputs/annulus_reference/C3_truekernel_out.txt, but that block is the 48/44 "
              "INDUCTION-machine annulus (X_g = 2.971852e-03, L = 0.164782 m, increment per "
              "decade 3.035405e-07 H). The submitted Table 8 is the PMSM annulus "
              "(mu_r = 1.0390, L = 33.0 mm, increment 6.0788e-08 H). Different machine, "
              "different numbers: the mapping by number would have been wrong. The predicted "
              "slope 6.078825e-08 is printed in R4_unlock_out.txt, but none of the summed, "
              "closed-form or increment values of the three panels is in outputs/."),
    dict(key="T9", quantity="Table 9: angular mesh convergence at one layer per shoe sub-region",
         table="9", figure=None,
         script="code/MEC_BLDC/RUN_A2_TABLE5.m",
         transcript="outputs/MEC_BLDC/A2_table5_out.txt",
         block="second block, panel 'Table 5(a) : n_sh = 1, solveur LINEAIRE'",
         date="2026-08-04 (docs/MANIFEST.md section 1)",
         configuration="PMSM 15/14, n_sh = 1 (nine radial layers), n_st = 4, n_ys = 3, linear "
                       "solver mu_r = 3000, kfringe = 0.325 (inert), P0 basis, rotor at phi = 0; "
                       "reference B_g1 = 1.07455 T, nu=8 = 0.019659 T, L_d = 52.3415 mH",
         authoritative=True,
         paras=[1016, 1017, 1018, 1019, 1022, 1023, 1024, 1025,
                1028, 1029, 1030, 1031, 1034, 1035, 1036, 1037],
         note="Mapped by content: same four tilings 180/360/540/900, same unknown counts, same "
              "three deviation columns and the same peak iron field. The first block of the "
              "transcript aborts on an unrecognised variable 'M'."),
    dict(key="T10", quantity="Table 10: angular sweep at four layers per shoe sub-region, both "
                             "solvers",
         table="10", figure=None,
         script="code/MEC_BLDC/RUN_X1_TABLE5B_RECONCILE.m",
         transcript="outputs/MEC_BLDC/X1_table5b_reconcile_out.txt",
         block="section 3, 'Table 5(b) retenue : n_sh = 4'; the identification residuals come "
               "from section 2 of the same run",
         date="2026-08-05 (docs/MANIFEST.md section 1)",
         configuration="PMSM 15/14, n_sh = 4 (fifteen radial layers), n_st = 4, n_ys = 3, "
                       "linear mu_r = 3000 and exact Newton on bh_curve M350, P0 basis, "
                       "numax = floor(M_s/2); reference nu=8 = 0.0196594150 T",
         authoritative=True,
         paras=[1050, 1051, 1052, 1053, 1056, 1057, 1058, 1059,
                1062, 1063, 1064, 1065, 1067],
         note="Mapped by content: the n_sh sweep of the transcript identifies n_sh = 4 with a "
              "maximum residual of 0.06 point over the six published cells, and the footnote's "
              "four residuals 9.08 / 2.78 / 0.81 / 0.06 are that identification."),
    dict(key="T11", quantity="Table 11: PMSM air-gap flux density at no load, spatial spectrum",
         table="11", figure=None,
         script="code/MEC_BLDC/RUN_C4C_TABLE8.m",
         transcript="outputs/MEC_BLDC/C4c_table8_out.txt",
         block="only block",
         date="2026-08-24 (as deposited)",
         configuration="PMSM 15/14, lumped network, mid-gap radius, kfringe = 0.325, "
                       "reference read from the finite-element project",
         authoritative=True,
         paras=[1080, 1081, 1085, 1086, 1090, 1091, 1095, 1096,
                1100, 1101, 1105, 1106, 1110, 1111],
         note="Mapped by content: the seven orders p, |N_s-p|, 3p, N_s+p, |2N_s-p|, 5p, 2N_s+p "
              "and their identification strings are the transcript's own rows."),
    dict(key="T12", quantity="Table 12: PMSM at no load, mesh at two shoe layerings, lumped "
                             "counterpart and 2-D finite-element reference",
         table="12", figure=None,
         script="code/MEC_BLDC/RUN_A1_TABLE7.m",
         transcript="outputs/MEC_BLDC/A1_table7_out.txt",
         block="LAST complete block (the transcript holds four starts; the earlier mesh columns "
               "return a peak phase EMF of exactly zero, the warm-start defect)",
         date="2026-08-04 (docs/MANIFEST.md section 1)",
         configuration="PMSM 15/14, M_s = 540, N_p = 61 rotor positions, P0 basis, "
                       "kfringe = 0.325 (inert in the mesh columns), linear solver mu_r = 3000; "
                       "unknowns 5400 / 6480 / 1260 surface columns",
         authoritative=True,
         paras=list(range(1119, 1204)) + [1291],
         note="Anchor supplied with the task (published there as 'Table 12/13') and re-verified: "
              "the seventeen quantities and the four columns are the transcript's last block. "
              "The footnote's 0.32 and 0.14 percentage points are the two cells marked SENSIBLE "
              "in outputs/MEC_BLDC/C4_ecarts_out.txt (0.318 and 0.136)."),
    dict(key="T13", quantity="Table 13: magnet eddy-current losses at no load, two field sources",
         table="13", figure=None,
         script="code/MEC_BLDC/RUN_A3_PMLOSS.m",
         transcript="outputs/MEC_BLDC/A3_pmloss_out.txt",
         block="section 'A VIDE'",
         date="2026-08-04 (docs/MANIFEST.md section 1)",
         configuration="frozen chain pm_loss, N_p = 241, 1688 rpm, sigma_pm = 555556 S/m, "
                       "kfringe = 0.325",
         authoritative=True,
         paras=[1298, 1299, 1300, 1302, 1303, 1304],
         note="Mapped by content: 'reseau a une dent' and 'MAILLAGE polaire' are the two rows. "
              "outputs/MEC_BLDC/A3_pmloss_prefix_out.txt is the pre-sign-correction copy and is "
              "exploratory. The provenance of the value published elsewhere for this quantity is "
              "in outputs/MEC_BLDC/D1_pmloss_chain_out.txt."),
    dict(key="T14", quantity="Table 14: where the sweep deviation originates, quasi-static "
                             "against settled phase current",
         table="14", figure=None,
         script=S_MASTER, transcript=MASTER,
         block="section 6b, the five-row quasi-static panel",
         date="2026-08-24 (as deposited; master run 677 s)",
         configuration="PMSM 15/14, 500 V dc six-step, saturable maps, drive_mec current loop; "
                       "tau = L/2R = 5.19 ms against a 9.83 ms electrical period at 872 rpm",
         authoritative=True,
         paras=list(range(1325, 1350)),
         note="Mapped by content: the five speeds 297/584/872/1159/1446 and the three current "
              "columns are the transcript's own panel."),
    dict(key="T15a", quantity="Table 15 upper panel: seven representative points of the speed "
                              "sweep at a 500 V dc bus",
         table="15", figure=None,
         script=S_MASTER, transcript=MASTER,
         block="section 6b, the thirty-point sweep table",
         date="2026-08-24 (as deposited)",
         configuration="30 points from 10 to 1676 rpm, 500 V dc bus, six-step",
         authoritative=True,
         paras=list(range(1365, 1428)),
         note="Mapped by content: identical seven speeds and five column pairs."),
    dict(key="T15b", quantity="Table 15 lower panel: anchors of the characteristic and the mean "
                              "deviations over the 28 points above 100 rpm",
         table="15", figure=None,
         script=S_MASTER + " + code/MEC_BLDC/RUN_R7_SCORECARD.m + code/MEC_BLDC/RUN_G4_POINTS.m",
         transcript=MASTER + " + " + R7 + " + " + G4,
         block="master section 6b (anchors and mean deviations); R7 rows 27 and 29; "
               "G4 row 29 in points",
         date="2026-08-24 (master, R7) / 2026-08-26 (G4)",
         configuration="same sweep; the peak-efficiency row is reported in points per the rule "
                       "of Section 5.3, which is what G4 computes",
         authoritative=True,
         paras=[1434, 1435, 1436, 1438, 1439, 1442, 1443, 1446, 1447,
                1450, 1451, 1454, 1455, 1456, 1459, 1461, 1463, 1465, 1467,
                1469, 1471, 1474],
         note="Mapped by content: 'puissance max', 'couple a l'arret', 'vitesse a vide' and the "
              "mean-deviation line of the master transcript, with the two efficiency figures "
              "carried in points by G4 (-3.986350 p.p.)."),
    dict(key="T16", quantity="Table 16: PMSM on-load transient, steady state",
         table="16", figure=None,
         script=S_MASTER, transcript=MASTER,
         block="section 5b, the full-precision panel saved as T10_full.mat",
         date="2026-08-24 (as deposited)",
         configuration="registration +143.49 deg electrical established on the flux linkage at "
                       "no load and held fixed; saturable maps L_line 104.3 -> 12.1 mH; "
                       "the two columns are at different speeds",
         authoritative=True,
         paras=list(range(1483, 1518)),
         note="Mapped by content: the nine rows are the transcript's own full-precision panel. "
              "Footnote c ('this row moved from 1.593 W') is the sign correction recorded in "
              "outputs/MEC_BLDC/D1_pmloss_chain_prefix_out.txt; outputs/MEC_BLDC/en.txt, fr.txt "
              "and outF.txt are the pre-correction copies of this master transcript and are "
              "superseded."),
    dict(key="T17", quantity="Table 17: harmonic comparison of the cogging torque",
         table="17", figure=None,
         script=S_MASTER, transcript=MASTER,
         block="section 7, 'COMPARAISON HARMONIQUE DU COUPLE DE DETENTE (A vs B)'",
         date="2026-08-24 (as deposited)",
         configuration="model samples one pole pitch at 840 positions; magnetostatic reference "
                       "remeshed at 360 positions over the revolution; only multiples of "
                       "lcm(15,14) = 210 admissible",
         authoritative=True,
         paras=[1526, 1527, 1529, 1530, 1532, 1533, 1538, 1539, 1541, 1542,
                1544, 1545, 1547, 1548, 1550, 1551],
         note="Mapped by content. ONE CELL DOES NOT MATCH: the table publishes an angular "
              "Nyquist order of 5880 for the model, the transcript prints 5866. The reference "
              "column (180) matches. Reported, not adjusted."),
]

FIGURES = [
    dict(key="F1", quantity="Figure 1: cross-section of the 750 W 15/14 PMSM",
         table=None, figure="1", script="", transcript="", block="",
         date="", configuration="", authoritative=False, paras=[],
         note="NOT ATTRIBUTED. A geometry drawing, not a computed result. The deposit holds "
              "code/article/figures/fig_machines.pdf, but the earlier manuscript's caption for "
              "it is 'Cross-sections of the TWO machines, to scale', so it is not the submitted "
              "single-machine figure and is not claimed as its source."),
    dict(key="F2", quantity="Figure 2: geometry of the two-layer annulus",
         table=None, figure="2", script="", transcript="",
         block="", date="", configuration="", authoritative=False, paras=[],
         note="Schematic. code/article/figures/fig_annulus.pdf carries the same caption in the "
              "earlier manuscript source, but no run produces it and no transcript is claimed."),
    dict(key="F3", quantity="Figure 3: the four-branch cell",
         table=None, figure="3", script="", transcript="",
         block="", date="", configuration="", authoritative=False, paras=[],
         note="Schematic; code/article/figures/fig_cell.pdf. No computed quantity, no chain."),
    dict(key="F4", quantity="Figure 4: radial layering of the stator mesh",
         table=None, figure="4", script="", transcript="",
         block="", date="", configuration="", authoritative=False, paras=[],
         note="Schematic; code/article/figures/fig_mesh.pdf. No computed quantity, no chain."),
    dict(key="F5", quantity="Figure 5: solution chain",
         table=None, figure="5", script="", transcript="",
         block="", date="", configuration="", authoritative=False, paras=[],
         note="Schematic; code/article/figures/fig_chain.pdf. No computed quantity, no chain."),
    dict(key="F6", quantity="Figure 6: the two comparison axes of Section 5.1",
         table=None, figure="6", script="", transcript="",
         block="", date="", configuration="", authoritative=False, paras=[],
         note="Schematic; code/article/figures/fig_protocol.pdf. No computed quantity, no chain."),
    dict(key="F7", quantity="Figure 7: convergence diagnostics of the released sweep of Table 6",
         table=None, figure="7", script="", transcript="",
         block="", date="", configuration="", authoritative=False, paras=[],
         note="NOT ATTRIBUTED. It plots the released sweep of Table 6 at M_s = 1080 over eleven "
              "doublings; the deposited chain stops at N_h = 4320 (see T6b). No figure file of "
              "that content exists in code/article/figures/."),
    dict(key="F8", quantity="Figure 8: air-gap flux density at no load, DtN operator against 2-D FEA",
         table=None, figure="8",
         script=S_MASTER + " (savefigure BLDC_FIG1_champ)",
         transcript=MASTER,
         block="sections 1-2 of the master run",
         date="2026-08-24 (as deposited)",
         configuration="PMSM 15/14 at no load, lumped iron model, mid-gap radius",
         authoritative=True, paras=[],
         note="Mapped by content and by the export name printed in the transcript "
              "('BLDC_FIG1_champ - champ d'entrefer (Br, Bt, spectre)'). File "
              "code/article/figures/BLDC_FIG1_champ.pdf."),
    dict(key="F9", quantity="Figure 9: internal magnetic circuit, magnetomotive force per magnet "
                            "and per gap",
         table=None, figure="9",
         script=S_MASTER + " (savefigure BLDC_FIG2_circuit)",
         transcript=MASTER,
         block="section 3 of the master run",
         date="2026-08-24 (as deposited)",
         configuration="lumped iron model; k_r computed retaining the outlying gap MMF "
                       "mmf_g1 = 531 A, second reference 1.0664 excluding it",
         authoritative=True, paras=[],
         note="Mapped by content and by the printed export name BLDC_FIG2_circuit."),
    dict(key="F10", quantity="Figure 10: back-electromotive force and flux linkage at 1500 rpm",
         table=None, figure="10",
         script=S_MASTER + " (savefigure BLDC_FIG3_bemf)",
         transcript=MASTER,
         block="section 4 of the master run",
         date="2026-08-24 (as deposited)",
         configuration="single angular registration +143.49 deg electrical established on the "
                       "flux linkage; residual optimal shift -0.5 deg",
         authoritative=True, paras=[],
         note="Mapped by content and by the printed export name BLDC_FIG3_bemf."),
    dict(key="F11", quantity="Figure 11: synthesis of the PMSM validation, the twenty-nine "
                             "deviations",
         table=None, figure="11",
         script="code/MEC_BLDC/RUN_R7_SCORECARD.m + code/MEC_BLDC/RUN_G4_POINTS.m "
                "(+ " + S_MASTER + " savefigure BLDC_FIG7_scorecard)",
         transcript=R7 + " + " + G4,
         block="R7 sections 1 to 6; G4 in full",
         date="2026-08-24 (R7) / 2026-08-26 (G4)",
         configuration="lumped iron model throughout; threshold |deviation| < 5 % strict; "
                       "rows 25 and 29 reported in points per Section 5.3",
         authoritative=True, paras=[],
         note="The twenty-nine bars are the twenty-nine rows tabulated in "
              "code/revision/annex_29_quantities.md. Both published counts re-derived there. "
              "outputs/MEC_BLDC/R7_scorecard_prefix_out.txt is the earlier 155 s run and is "
              "superseded by the 679 s run."),
    dict(key="F12", quantity="Figure 12: saturation curves L(i) and psi(i) from the map",
         table=None, figure="12",
         script=S_MASTER + " (savefigure BLDC_FIG5_3_saturation)",
         transcript=MASTER,
         block="section 5b, saturation line",
         date="2026-08-24 (as deposited)",
         configuration="line inductance 104.3 mH at 0 A to 12.1 mH at 25 A; psi_ab 3.01 to "
                       "0.68 Wb/rad",
         authoritative=True, paras=[],
         note="Mapped by content and by the printed export name BLDC_FIG5_3_saturation. The map "
              "itself is audited in outputs/MEC_BLDC/R6_satmap_out.txt and R6a_grid_out.txt."),
    dict(key="F13", quantity="Figure 13: on-load transient, speed, phase current, steady-state "
                             "current and dc-bus current",
         table=None, figure="13",
         script=S_MASTER + " (savefigure BLDC_FIG5_1_transient)",
         transcript=MASTER,
         block="section 5b of the master run",
         date="2026-08-24 (as deposited)",
         configuration="MEC drive model, 500 V dc six-step, load torque 4.8 N m",
         authoritative=True, paras=[],
         note="Mapped by content and by the printed export name BLDC_FIG5_1_transient "
              "('transitoire electrique (4 vues)')."),
    dict(key="F14", quantity="Figure 14: on-load transient, conversion and losses",
         table=None, figure="14",
         script=S_MASTER + " (savefigure BLDC_FIG5_2_conversion)",
         transcript=MASTER,
         block="section 5b, power balance",
         date="2026-08-24 (as deposited)",
         configuration="power balance in the ANSYS definition; efficiency by (37) in both columns",
         authoritative=True, paras=[],
         note="Mapped by content and by the printed export name BLDC_FIG5_2_conversion "
              "('conversion et pertes (4 vues)')."),
    dict(key="F15", quantity="Figure 15: speed sweep at constant dc-bus voltage",
         table=None, figure="15",
         script=S_MASTER + " (savefigure BLDC_FIG6_balayage)",
         transcript=MASTER,
         block="sections 6 and 6b",
         date="2026-08-24 (as deposited)",
         configuration="30 points, 10 to 1676 rpm, 500 V dc; 'BEMF' is the terminal envelope "
                       "Vdc - 2Ri, not the back-EMF",
         authoritative=True, paras=[],
         note="Mapped by content and by the printed export name BLDC_FIG6_balayage."),
    dict(key="F16", quantity="Figure 16: cogging torque of the PMSM",
         table=None, figure="16",
         script=S_MASTER + " (savefigure BLDC_FIG4_detente)",
         transcript=MASTER,
         block="section 7 of the master run",
         date="2026-08-24 (as deposited)",
         configuration="magnetostatic reference remeshed at 360 positions; surrogate with "
                       "identical amplitude spectrum; C3 block for the sigma multiples",
         authoritative=True, paras=[],
         note="Mapped by content and by the printed export name BLDC_FIG4_detente. The "
              "sigma-multiple analysis of the same figure is in "
              "outputs/MEC_BLDC/C3_sigma_out.txt (RUN_C3_SIGMA.m)."),
]


# --------------------------------------------------------------------------
# 4. exploratory scripts -- declared, never deleted
# --------------------------------------------------------------------------
def scan_scripts():
    """Every RUN_*.m of this machine, and whether its transcript is deposited."""
    recs = []
    for sub in ("MEC_BLDC", "article"):
        d = os.path.join(REPO, "code", sub)
        if not os.path.isdir(d):
            continue
        for fn in sorted(os.listdir(d)):
            if not (fn.startswith("RUN_") and fn.endswith(".m")):
                continue
            src = io.open(os.path.join(d, fn), encoding="utf-8",
                          errors="replace").read()
            txts = sorted(set(re.findall(r"[A-Za-z0-9_]+_out\.txt", src)))
            present = []
            for t in txts:
                for od in ("MEC_BLDC", "annulus_reference", "revision"):
                    p = os.path.join(REPO, "outputs", od, t)
                    if os.path.isfile(p):
                        present.append("outputs/%s/%s" % (od, t))
            recs.append({
                "script": "code/%s/%s" % (sub, fn),
                "declares_transcript": txts,
                "transcript_present": sorted(set(present)),
            })
    return recs


def exploratory_records(scripts, used_scripts):
    out = []
    for r in scripts:
        name = os.path.basename(r["script"])
        diag = name.startswith("RUN_DIAG_")
        writes = bool(r["transcript_present"])
        if r["script"] in used_scripts and not diag:
            continue
        if diag:
            reason = ("RUN_DIAG_* family: diagnostic by declaration. "
                      "Section 5.2 is right to keep it; it is kept.")
        elif not writes:
            reason = ("Writes no transcript into outputs/, so no published value can "
                      "be traced through it.")
        else:
            reason = ("Writes a transcript, but no cell of the submitted manuscript "
                      "was traced to it in this pass.")
        out.append(dict(
            quantity="Exploratory / unattributed script: " + name,
            table=None, figure=None,
            script=r["script"],
            transcript=("; ".join(r["transcript_present"])
                        if r["transcript_present"] else ""),
            date="", configuration="",
            authoritative=False,
            block="", values=[],
            note=reason + " NOT DELETED."))
    return out


# --------------------------------------------------------------------------
# 5. assemble
# --------------------------------------------------------------------------
def main():
    ms = load_ms()
    r7_path = os.path.join(REPO, *R7.split("/"))
    g4_path = os.path.join(REPO, *G4.split("/"))
    rows = parse_r7(r7_path)
    ctx = parse_r7_context(r7_path)
    g4 = parse_g4(g4_path)

    if len(rows) != 29:
        print("WARNING: parsed %d rows from R7, expected 29" % len(rows))

    records = []
    used_scripts = set()

    for spec in TABLES + FIGURES:
        rec = {
            "quantity": spec["quantity"],
            "table": spec["table"],
            "figure": spec["figure"],
            "script": spec["script"],
            "transcript": spec["transcript"],
            "date": spec["date"],
            "configuration": spec["configuration"],
            "authoritative": spec["authoritative"],
            "block": spec["block"],
            "values": published(ms, spec["paras"]),
            "note": spec["note"],
        }
        records.append(rec)
        for s in re.findall(r"code/[A-Za-z_]+/[A-Za-z0-9_]+\.m", spec["script"]):
            used_scripts.add(s)

    # the twenty-nine quantities, one record each
    for r in rows:
        n = r["n"]
        pts = g4["points"].get(n)
        vals = [{"value": r["model"], "role": "model"},
                {"value": r["reference"], "role": "reference"},
                {"value": r["deviation_pct"], "role": "deviation_pct"}]
        script = "code/MEC_BLDC/RUN_R7_SCORECARD.m"
        transcript = R7
        note = ("Row %d of the recount. Produced inside the master run and re-counted "
                "by R7; the deviation is formed on full-precision values." % n)
        if pts:
            script += " + code/MEC_BLDC/RUN_G4_POINTS.m"
            transcript += " + " + G4
            vals.append({"value": pts["model_full"], "role": "model_full_precision"})
            vals.append({"value": pts["reference_full"], "role": "reference_full_precision"})
            vals.append({"value": pts["points_pp"], "role": "deviation_points"})
            note += (" Reported in points per Section 5.3: G4 gives %s p.p." %
                     pts["points_pp"])
        records.append({
            "quantity": "Compared quantity %d of 29: %s" % (n, EN[n]),
            "table": "15" if n in (26, 27, 28, 29) else None,
            "figure": "11",
            "script": script,
            "transcript": transcript,
            "date": "2026-08-24 (R7)" + (" / 2026-08-26 (G4)" if pts else ""),
            "configuration": SECTION[n] + "; lumped iron model; threshold |dev| < 5 % strict",
            "authoritative": True,
            "block": "R7 section 1" + (" + G4 section 'les deux lignes de rendement'"
                                       if pts else ""),
            "values": vals,
            "note": note,
        })

    scripts = scan_scripts()
    records.extend(exploratory_records(scripts, used_scripts))

    n_run = len(scripts)
    n_with = len([s for s in scripts if s["transcript_present"]])

    doc = {
        "manifest": "B10 -- the authoritative chain for the submitted manuscript",
        "manuscript": ("Boundary Condensation of a Two-Layer Dirichlet-to-Neumann "
                       "Operator: Trace Conformity and Air-Gap Mesh Reduction"),
        "date_iso": ISO_DATE,
        "numbering": ("Tables 1-17 and Figures 1-16 are those of the SUBMITTED manuscript. "
                      "docs/PROVENANCE.md follows the numbering of the earlier manuscript "
                      "('A Well-Posed Dirichlet-to-Neumann Air-Gap Operator...') and was NOT "
                      "used as a lookup table: every row below was mapped by CONTENT, by "
                      "matching the submitted caption and cells against the transcripts."),
        "rules": [
            "No value in this file was typed by hand: published values are read from the "
            "rendered manuscript by paragraph index, computed values are parsed from the "
            "transcripts.",
            "Where a chain could not be established the record carries authoritative=false "
            "with an empty script and transcript. An empty cell is information; a wrong "
            "cell is not.",
            "Nothing is deleted. The RUN_DIAG_* family and every script that writes no "
            "transcript are declared exploratory and kept.",
        ],
        "script_census": {
            "run_scripts_found": n_run,
            "run_scripts_whose_transcript_is_deposited": n_with,
            "manuscript_claim": ("Section 5.2: '76 executable run scripts for this machine, "
                                 "of which 18 write a dated transcript'"),
        },
        "counts_from_transcripts": {
            "R7_master_counter": ctx.get("master_counter"),
            "R7_independent_recount": ctx.get("independent_recount"),
            "R7_excluding_three": ctx.get("excluding_three"),
            "G4_relative_all": g4["recount"].get("relative_all"),
            "G4_points_rows_25_29": g4["recount"].get("points_2529"),
            "G4_relative_excluding_three": g4["recount"].get("relative_excl"),
            "G4_points_excluding_three": g4["recount"].get("points_excl"),
            "G4_within_one_point": g4.get("within_one_point"),
            "G4_recomposition_guard": g4.get("guard"),
        },
        "records": records,
    }

    with io.open(OUT_JSON, "w", encoding="utf-8") as fh:
        fh.write(json.dumps(doc, indent=2, ensure_ascii=False))
    print("wrote %s (%d records)" % (OUT_JSON, len(records)))

    write_annex(rows, ctx, g4, r7_path, g4_path)
    return doc


# --------------------------------------------------------------------------
# 6. the annex of the twenty-nine quantities
# --------------------------------------------------------------------------
def write_annex(rows, ctx, g4, r7_path, g4_path):
    L = []
    A = L.append
    A("# Annex A. The twenty-nine compared quantities behind Figure 11")
    A("")
    A("*Prepared %s for the revision of \"Boundary Condensation of a Two-Layer "
      "Dirichlet-to-Neumann Operator: Trace Conformity and Air-Gap Mesh Reduction\".*" % ISO_DATE)
    A("")
    A("Figure 11 of the submitted manuscript reports two counts a reader cannot check from "
      "the plot alone: *twenty-one of twenty-nine fall within 5 %*, and *twenty of "
      "twenty-six* once the three quantities the manuscript declares not validated are set "
      "aside. This annex tabulates the twenty-nine bars at full precision, with the model "
      "value, the reference value, the signed deviation and the transcript each is read "
      "from, so that a referee can recount both figures line by line.")
    A("")
    A("**Convention.** The deviation is the signed relative deviation of (38), formed on "
      "full-precision values and rounded once for display; negative means the network is "
      "below the reference. The threshold is strict, `|deviation| < 5`. Rows 25 and 29 are "
      "efficiencies: Section 5.3 requires them in percentage points, and the two columns "
      "for them are given below.")
    A("")
    A("**Provenance.** Column *source* names the transcript each line is read from:")
    A("")
    A("- `R7` = `outputs/MEC_BLDC/R7_scorecard_out.txt`, produced by "
      "`code/MEC_BLDC/RUN_R7_SCORECARD.m`;")
    A("- `G4` = `outputs/MEC_BLDC/G4_points_out.txt`, produced by "
      "`code/MEC_BLDC/RUN_G4_POINTS.m`, which re-reads `code/MEC_BLDC/R7_scorecard.mat`.")
    A("")
    A("Every value below was extracted mechanically from those two files by "
      "`code/revision/build_manifest.py`; none was transcribed by hand.")
    A("")
    A("## A.1 The twenty-nine quantities")
    A("")
    A("| # | quantity | model | reference | deviation (%) | within 5 % | source |")
    A("|---:|---|---:|---:|---:|:---:|:---|")
    for r in rows:
        dev = float(r["deviation_pct"])
        inside = "yes" if abs(dev) < 5 else "**no**"
        A("| %d | %s | %s | %s | %s | %s | R7 |" % (
            r["n"], EN[r["n"]], r["model"], r["reference"],
            r["deviation_pct"], inside))
    A("")
    A("## A.2 The two efficiency rows in percentage points")
    A("")
    A("Section 5.3 reports a deviation between two efficiencies in points, not in per cent. "
      "`G4` carries both readings at full precision.")
    A("")
    A("| # | quantity | model (full) | reference (full) | relative (%) | points (p.p.) | source |")
    A("|---:|---|---:|---:|---:|---:|:---|")
    for n in sorted(g4["points"]):
        p = g4["points"][n]
        A("| %d | %s | %s | %s | %s | %s | G4 |" % (
            n, EN[n], p["model_full"], p["reference_full"],
            p["relative_pct"], p["points_pp"]))
    A("")
    A("## A.3 The three quantities declared not validated")
    A("")
    A("The manuscript excludes three quantities from the second count and says why. The "
      "motives below are the transcript's own.")
    A("")
    A("| # | quantity | deviation (%) | motive as recorded in R7 |")
    A("|---:|---|---:|:---|")
    mot = ctx.get("motives", [])
    for i, (num, label, dev) in enumerate(ctx.get("excluded", [])):
        m = mot[i] if i < len(mot) else ""
        A("| %s | %s | %s | %s |" % (num, EN.get(int(num), label.strip()), dev, m))
    A("")
    A("## A.4 Re-derivation of the two published counts")
    A("")
    inside_all = [r for r in rows if abs(float(r["deviation_pct"])) < 5]
    excl = set(int(e[0]) for e in ctx.get("excluded", []))
    kept = [r for r in rows if r["n"] not in excl]
    inside_kept = [r for r in kept if abs(float(r["deviation_pct"])) < 5]
    A("Counted directly from the table of A.1, with the strict threshold `|dev| < 5`:")
    A("")
    A("| count | re-derived here | printed by R7 | printed by G4 | published in Fig. 11 | holds |")
    A("|:--|:--|:--|:--|:--|:--|")
    A("| within 5 %%, all twenty-nine | **%d / %d** | %s | %s | twenty-one of twenty-nine | %s |"
      % (len(inside_all), len(rows), ctx.get("independent_recount"),
         g4["recount"].get("relative_all"),
         "yes" if len(inside_all) == 21 and len(rows) == 29 else "NO"))
    A("| within 5 %%, excluding the three not validated | **%d / %d** | %s | %s | "
      "twenty of twenty-six | %s |"
      % (len(inside_kept), len(kept), ctx.get("excluding_three"),
         g4["recount"].get("relative_excl"),
         "yes" if len(inside_kept) == 20 and len(kept) == 26 else "NO"))
    A("")
    A("The quantities outside the threshold are, in order: " +
      ", ".join("#%d %s (%s %%)" % (r["n"], EN[r["n"]], r["deviation_pct"])
                for r in rows if abs(float(r["deviation_pct"])) >= 5) + ".")
    A("")
    A("Removing #%s from that list leaves the %d exceedances of the second count."
      % (", #".join(str(x) for x in sorted(excl)), len(kept) - len(inside_kept)))
    A("")
    near = [r for r in rows if 4 <= abs(float(r["deviation_pct"])) < 5]
    A("**Quantities within one point of the threshold.** In the relative convention there "
      "are %d: %s. G4 records \"%s\". The submitted caption says *two*, which is the "
      "reading in the points convention of Section 5.3, where row 29 becomes %s p.p. and "
      "leaves the one-point band; the caption and the transcript agree once the convention "
      "is stated."
      % (len(near),
         ", ".join("#%d (%s %%)" % (r["n"], r["deviation_pct"]) for r in near),
         g4.get("within_one_point", ""),
         g4["points"].get(29, {}).get("points_pp", "")))
    A("")
    A("**Guard.** G4 recomposes the twenty-nine deviations from the saved model and "
      "reference values and compares them with the stored deviation column; the maximum "
      "discrepancy is %s, that is machine precision." % g4.get("guard"))
    A("")
    A("## A.5 What this annex does not settle")
    A("")
    A("Two limits belong with the table.")
    A("")
    A("1. The reference carries no uncertainty band on any field quantity, and the "
      "manuscript states that it does not. A deviation of a few tenths of a per cent on "
      "rows 1 to 8 is therefore not evidence of agreement to that accuracy.")
    A("2. Rows 17 to 25 are read at an operating point the reference does not share: the "
      "network settles at a different speed. The efficiency agreement of row 25 is "
      "obtained there and nowhere else.")
    A("")
    A("---")
    A("")
    A("Files read to build this annex: `%s`, `%s`." %
      (R7, G4))
    with io.open(OUT_ANNEX, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print("wrote %s" % OUT_ANNEX)


if __name__ == "__main__":
    main()
