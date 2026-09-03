# -*- coding: utf-8 -*-
"""B12 - terminology and status consistency checker for the submitted manuscript.

Runs the ten checks C-01 .. C-10 directly on word/document.xml inside the .docx
(zipfile + ElementTree only; no python-docx, no scipy).

Purpose is structural, not cosmetic: it exists so that a partially corrected
manuscript FAILS.  A check that cannot be settled mechanically reports
UNDECIDED and prints its candidate sentences; it never reports a silent PASS.

Exit code 0 only if every check is PASS.  Any FAIL or UNDECIDED gives 1.

Usage:
    python consistency_check.py [path\\to\\manuscript.docx] [--out FILE]

Pure standard library + numpy.
"""

from __future__ import annotations

import argparse
import io
import os
import re
import sys
import zipfile
import unicodedata
import xml.etree.ElementTree as ET
from collections import OrderedDict, defaultdict

import numpy as np  # used for the sequence / coverage arithmetic

# --------------------------------------------------------------------------
# constants of this task
# --------------------------------------------------------------------------

ISO_DATE = "2026-09-02"

DEFAULT_DOCX = (r"C:\Users\hp\Desktop\Boundary Condensation of a Two-Layer "
                r"Dirichlet-to-Neumann Operator Trace Conformity and "
                r"Air-Gap Mesh Reduction.docx")

DEFAULT_OUT = (r"C:\Users\hp\Desktop\mec-dtn-airgap-pmsm\outputs\revision"
               r"\B12_consistency_check_out.txt")

W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
M = '{http://schemas.openxmlformats.org/officeDocument/2006/math}'

# --------------------------------------------------------------------------
# 1.  document model
# --------------------------------------------------------------------------

# Unicode sub/superscript letters and digits, mapped back to ASCII so that
# "Nₕ", "Mₛ", "Δₖ₊₁" and "10⁻⁸" become "N_{h}", "M_{s}", "Delta_{k+1}", "10^{-8}".
_SUBS = {
    '\u2080': '0', '\u2081': '1', '\u2082': '2', '\u2083': '3', '\u2084': '4',
    '\u2085': '5', '\u2086': '6', '\u2087': '7', '\u2088': '8', '\u2089': '9',
    '\u208a': '+', '\u208b': '-', '\u208c': '=', '\u208d': '(', '\u208e': ')',
    '\u2090': 'a', '\u2091': 'e', '\u2092': 'o', '\u2093': 'x', '\u2095': 'h',
    '\u2096': 'k', '\u2097': 'l', '\u2098': 'm', '\u2099': 'n', '\u209a': 'p',
    '\u209b': 's', '\u209c': 't', '\u1d62': 'i', '\u1d63': 'r', '\u1d64': 'u',
    '\u1d65': 'v', '\u2c7c': 'j',
}
_SUPS = {
    '\u2070': '0', '\u00b9': '1', '\u00b2': '2', '\u00b3': '3', '\u2074': '4',
    '\u2075': '5', '\u2076': '6', '\u2077': '7', '\u2078': '8', '\u2079': '9',
    '\u207a': '+', '\u207b': '-', '\u207c': '=', '\u207d': '(', '\u207e': ')',
    '\u207f': 'n', '\u1d43': 'a', '\u1d47': 'b', '\u1d9c': 'c', '\u1d48': 'd',
    '\u1d49': 'e', '\u1da0': 'f', '\u1d4d': 'g', '\u02b0': 'h', '\u2071': 'i',
    '\u02b2': 'j', '\u1d4f': 'k', '\u02e1': 'l', '\u1d50': 'm', '\u1d52': 'o',
    '\u1d56': 'p', '\u02b3': 'r', '\u02e2': 's', '\u1d57': 't', '\u1d58': 'u',
    '\u1d5b': 'v', '\u02b7': 'w', '\u02e3': 'x', '\u02b8': 'y', '\u1d41': 'A',
    '\u1d2c': 'A', '\u1d2e': 'B', '\u1d30': 'D', '\u1d31': 'E', '\u1d33': 'G',
    '\u1d34': 'H', '\u1d35': 'I', '\u1d36': 'J', '\u1d37': 'K', '\u1d38': 'L',
    '\u1d39': 'M', '\u1d3a': 'N', '\u1d3c': 'O', '\u1d3e': 'P', '\u1d3f': 'R',
    '\u1d40': 'T', '\u1d41': 'U', '\u1d42': 'W', '\u1d5a': 'U',
}
_SUB_RE = re.compile('[' + ''.join(_SUBS) + ']+')
_SUP_RE = re.compile('[' + ''.join(_SUPS) + ']+')


def normalise_scripts(s: str) -> str:
    """Turn runs of Unicode sub/superscript characters into _{...} / ^{...}."""
    if not s:
        return s
    s = _SUB_RE.sub(lambda m: '_{' + ''.join(_SUBS[c] for c in m.group(0)) + '}', s)
    s = _SUP_RE.sub(lambda m: '^{' + ''.join(_SUPS[c] for c in m.group(0)) + '}', s)
    return s


def _omml(el) -> str:
    """Linearise one OMML subtree, keeping sub/superscript structure in braces."""
    t = el.tag
    if t in (M + 'r', W + 'r'):
        return ''.join((x.text or '') for x in el.iter() if x.tag in (M + 't', W + 't'))
    if t == M + 'sSub':
        b = ''.join(_omml(c) for c in el if c.tag == M + 'e')
        s = ''.join(_omml(c) for c in el if c.tag == M + 'sub')
        return '%s_{%s}' % (b, s)
    if t == M + 'sSup':
        b = ''.join(_omml(c) for c in el if c.tag == M + 'e')
        s = ''.join(_omml(c) for c in el if c.tag == M + 'sup')
        return '%s^{%s}' % (b, s)
    if t == M + 'sSubSup':
        b = ''.join(_omml(c) for c in el if c.tag == M + 'e')
        sb = ''.join(_omml(c) for c in el if c.tag == M + 'sub')
        sp = ''.join(_omml(c) for c in el if c.tag == M + 'sup')
        return '%s_{%s}^{%s}' % (b, sb, sp)
    if t == M + 'f':
        n = ''.join(_omml(c) for c in el if c.tag == M + 'num')
        d = ''.join(_omml(c) for c in el if c.tag == M + 'den')
        return '(%s)/(%s)' % (n, d)
    if t == M + 'd':
        return '(%s)' % ''.join(_omml(c) for c in el if c.tag == M + 'e')
    if t == M + 'rad':
        return 'sqrt(%s)' % ''.join(_omml(c) for c in el if c.tag == M + 'e')
    if t == M + 'nary':
        sb = ''.join(_omml(c) for c in el if c.tag == M + 'sub')
        sp = ''.join(_omml(c) for c in el if c.tag == M + 'sup')
        e = ''.join(_omml(c) for c in el if c.tag == M + 'e')
        return '\u2211_{%s}^{%s}[%s]' % (sb, sp, e)
    out = []
    for c in el:
        if c.tag.endswith('Pr'):
            continue
        out.append(_omml(c))
    return ''.join(out)


class Para(object):
    __slots__ = ('idx', 'style', 'in_table', 'text', 'maths', 'numid',
                 'section', 'bucket', 'is_caption', 'is_note')

    def __init__(self, idx, style, in_table, text, maths, numid):
        self.idx = idx
        self.style = style
        self.in_table = in_table
        self.text = text
        self.maths = maths
        self.numid = numid
        self.section = ''
        self.bucket = ''
        self.is_caption = False
        self.is_note = False

    def __repr__(self):
        return '<p%d %s %s>' % (self.idx, self.style, self.text[:40])


class Manuscript(object):
    def __init__(self, path):
        self.path = path
        z = zipfile.ZipFile(path)
        self.root = ET.fromstring(z.read('word/document.xml'))
        try:
            self.numbering = ET.fromstring(z.read('word/numbering.xml'))
        except KeyError:
            self.numbering = None
        z.close()
        self.body = self.root.find(W + 'body')
        self._read_paragraphs()
        self._number_sections()
        self._classify()

    # -- paragraphs -------------------------------------------------------
    def _read_paragraphs(self):
        tbl_ids = set()
        self.tables = []
        for t in self.body.iter(W + 'tbl'):
            rows = []
            for tr in t.findall(W + 'tr'):
                cells = []
                for tc in tr.findall(W + 'tc'):
                    ps = list(tc.iter(W + 'p'))
                    cells.append(ps)
                rows.append(cells)
            self.tables.append(rows)
            for p in t.iter(W + 'p'):
                tbl_ids.add(id(p))

        self.paras = []
        self._by_id = {}
        for i, p in enumerate(self.body.iter(W + 'p'), 1):
            oms = list(p.iter(M + 'oMath'))
            inside = set()
            for om in oms:
                for d in om.iter():
                    inside.add(id(d))
            pieces, maths = [], []
            for node in p.iter():
                if node.tag == M + 'oMath':
                    s = normalise_scripts(_omml(node))
                    maths.append(s)
                    pieces.append(' ' + s + ' ')
                elif node.tag == W + 't' and id(node) not in inside:
                    pieces.append(normalise_scripts(node.text or ''))
            st = p.find(W + 'pPr/' + W + 'pStyle')
            style = st.get(W + 'val') if st is not None else ''
            nm = p.find(W + 'pPr/' + W + 'numPr/' + W + 'numId')
            numid = nm.get(W + 'val') if nm is not None else None
            text = re.sub(r'[ \t\u00a0]+', ' ', ''.join(pieces)).strip()
            par = Para(i, style or 'body', id(p) in tbl_ids, text, maths, numid)
            self.paras.append(par)
            self._by_id[id(p)] = par

    def cell_paras(self, cell):
        return [self._by_id[id(p)] for p in cell if id(p) in self._by_id]

    # -- section numbering ------------------------------------------------
    BACK_MATTER = ('statements and declarations', 'data availability',
                   'references', 'acknowledgements', 'acknowledgments',
                   'appendix', 'nomenclature', 'abbreviations')

    def _number_sections(self):
        h1 = h2 = h3 = 0
        cur = ''
        started = False
        for p in self.paras:
            if p.style == 'Heading1':
                if p.text.strip().lower() in self.BACK_MATTER:
                    cur = 'back'
                    h2 = h3 = 0
                else:
                    h1 += 1
                    h2 = h3 = 0
                    cur = str(h1)
                    started = True
            elif p.style == 'Heading2' and cur not in ('', 'back'):
                h2 += 1
                h3 = 0
                cur = '%d.%d' % (h1, h2)
            elif p.style == 'Heading3' and cur not in ('', 'back') and h2:
                h3 += 1
                cur = '%d.%d.%d' % (h1, h2, h3)
            p.section = cur if started or cur == 'back' else 'front'
        # the set of section labels that exist
        self.sections = OrderedDict()
        for p in self.paras:
            if p.style in ('Heading1', 'Heading2', 'Heading3') and \
                    p.section not in ('back', 'front'):
                self.sections.setdefault(p.section, p.text.strip())

    CAPTION_RE = re.compile(r'^\s*(Table|Figure|Fig\.)\s+(\d+)\s*[.:]')

    def _classify(self):
        # abstract block: between the paragraph 'Abstract' and 'Keywords'
        self.abstract_paras = []
        seen_abs = False
        for p in self.paras[:40]:
            low = p.text.strip().lower()
            if low == 'abstract' or low.startswith('abstract '):
                seen_abs = True
                if low != 'abstract':
                    self.abstract_paras.append(p)
                continue
            if seen_abs:
                if low.startswith('keywords') or low.startswith('key words'):
                    break
                if p.text.strip():
                    self.abstract_paras.append(p)
        abs_ids = set(id(p) for p in self.abstract_paras)

        for p in self.paras:
            m = self.CAPTION_RE.match(p.text)
            p.is_caption = bool(m)
            # a "table note" is caption prose or a long prose paragraph living
            # inside a w:tbl (footnote rows of the big tables)
            p.is_note = bool(m) or (p.in_table and len(p.text) > 140)
            if id(p) in abs_ids:
                p.bucket = 'Abstract'
            elif p.is_note:
                p.bucket = 'table notes'
            elif p.section in ('back', 'front', ''):
                p.bucket = ''
            else:
                p.bucket = 'S' + p.section


# --------------------------------------------------------------------------
# 2.  sentence splitting
# --------------------------------------------------------------------------

_ABBR = ['Fig', 'Figs', 'Eq', 'Eqs', 'Sec', 'Sect', 'Tab', 'Ref', 'Refs', 'No',
         'vs', 'cf', 'approx', 'e.g', 'i.e', 'Dr', 'Prof', 'St', 'al', 'ca',
         'Univ', 'Press', 'ed', 'vol', 'pp', 'Inst', 'Trans', 'Mr', 'Ms']
_ABBR_RE = re.compile(r'\b(' + '|'.join(re.escape(a) for a in _ABBR) + r')\.')
_SENT_SPLIT = re.compile(r'(?<=[.!?])[ \u00a0]+(?=[A-Z\u00c0-\u00de(\u201c"\u2018])')


def sentences(text):
    if not text.strip():
        return []
    guard = _ABBR_RE.sub(lambda m: m.group(1) + '\x00', text)
    # protect decimals and numbered items such as "1.0390" or "Table 12."
    guard = re.sub(r'(?<=\d)\.(?=\d)', '\x00', guard)
    parts = _SENT_SPLIT.split(guard)
    return [p.replace('\x00', '.').strip() for p in parts if p.strip()]


def clip(s, n=300):
    s = ' '.join(s.split())
    return s if len(s) <= n else s[:n - 3] + '...'


# --------------------------------------------------------------------------
# 3.  reporting scaffolding
# --------------------------------------------------------------------------

class Report(object):
    def __init__(self):
        self.buf = io.StringIO()
        self.verdicts = OrderedDict()

    def w(self, s=''):
        self.buf.write(s + '\n')

    def verdict(self, cid, title, result, detail=''):
        self.verdicts[cid] = (title, result, detail)

    def text(self):
        return self.buf.getvalue()


# --------------------------------------------------------------------------
# 4.  the ten checks
# --------------------------------------------------------------------------

# ---- C-01 -----------------------------------------------------------------

SWEEP = ['Abstract', 'S1', 'S6.5', 'S7.2', 'S7.3', 'S8', 'table notes']

QUANTITIES = OrderedDict([
    ('peak flux linkage',
     r'(peak\s+)?flux linkage'),
    ('L_d (synchronous inductance)',
     r'synchronous inductance|\bL_\{d\}|\bLd\b'),
    ('L_a (phase self-inductance)',
     r'(phase\s+)?self-?inductance|\bL_\{a\}'),
    ('M (mutual inductance)',
     r'mutual inductance|\bthe mutual\b'),
    ('B_g1 (fundamental gap flux density)',
     r'fundamental (air-)?gap flux density|fundamental of the gap flux density|'
     r'fundamental air-gap flux density|B_\{g1\}'),
    ('nu = 8 sideband',
     r'\u03bd\s*=\s*8|first slot sideband|first sideband|slot sidebands?'),
    ('cogging torque',
     r'cogging'),
    ('magnet loss',
     r'magnet (eddy-current )?loss|magnet losses'),
    ('standstill torque',
     r'standstill torque'),
])

STATUS_CUES = OrderedDict([
    ('not validated',
     r'\bnot validated\b|\bdeclared not validated\b|\bcannot be validated\b|'
     r'\bis not a validation\b|\bnon-rejection\b'),
    ('accuracy',
     r'\baccuracy\b|\baccurate(?:ly)?\b'),
    ('validated',
     r'\bvalidat(?:ed|ion|es|ing)\b|\bagreement\b|\bagrees\b'),
    ('bound',
     r'\bbounds?\b|\bbounded\b|\bbracket(?:s|ed|ing)?\b|\bfloor\b|'
     r'\bno larger than the reference can resolve\b'),
    ('consistency indicator',
     r'consistency indicator|internal consistency|consistency check|'
     r'\bdiagnostic\b'),
])

NEGATORS = (r'\bnot\b|\bno\b|\bnone\b|\bnever\b|\bneither\b|\bnothing\b|'
            r'\bcannot\b|\bcan not\b|\brather than\b|\binstead of\b|'
            r'\bwithout\b|\bfails? to\b|\bis not\b|\bare not\b')
_NEG_RE = re.compile(NEGATORS, re.I)

EXCLUSIVE = {'accuracy', 'validated', 'bound', 'not validated',
             'consistency indicator'}


def status_of(sentence):
    """Return (positive_statuses, negated_statuses) found in one sentence."""
    pos, neg = set(), set()
    for name, pat in STATUS_CUES.items():
        for m in re.finditer(pat, sentence, re.I):
            left = sentence[max(0, m.start() - 60):m.start()]
            if _NEG_RE.search(left):
                neg.add(name)
            else:
                pos.add(name)
    # "not validated" subsumes the "validated" token it contains
    if 'not validated' in pos:
        pos.discard('validated')
        neg.add('validated')
    return pos, neg


def check_c01(ms, rep):
    rep.w('-' * 78)
    rep.w('[C-01]  status of every compared quantity, section by section')
    rep.w('-' * 78)
    rep.w('Sweep: ' + ', '.join(SWEEP))
    rep.w('Status vocabulary: accuracy | validated | bound | consistency '
          'indicator | not validated')
    rep.w('A cue preceded by a negator within 60 characters is recorded as '
          'DENIED, not as a status.')
    rep.w()

    table = defaultdict(lambda: defaultdict(lambda: [set(), set(), []]))
    for p in ms.paras:
        if p.bucket not in SWEEP:
            continue
        for sent in sentences(p.text):
            for qname, qpat in QUANTITIES.items():
                if not re.search(qpat, sent, re.I):
                    continue
                pos, neg = status_of(sent)
                cell = table[qname][p.bucket]
                cell[0] |= pos
                cell[1] |= neg
                if pos or neg:
                    cell[2].append((p.idx, sent))
                elif len(cell[2]) < 12:
                    cell[2].append((p.idx, sent))

    conflicts, silent = [], []
    for qname in QUANTITIES:
        rep.w('QUANTITY: %s' % qname)
        seen_pos = set()
        any_section = False
        for bucket in SWEEP:
            if bucket not in table[qname]:
                rep.w('    %-13s  --  (not mentioned)' % bucket)
                continue
            any_section = True
            pos, neg, sents = table[qname][bucket]
            seen_pos |= (pos & EXCLUSIVE)
            lab = ', '.join(sorted(pos)) if pos else '(no status word)'
            if neg:
                lab += '   [DENIED: %s]' % ', '.join(sorted(neg))
            rep.w('    %-13s  ->  %s' % (bucket, lab))
            for idx, s in sents[:4]:
                rep.w('         p%-5d "%s"' % (idx, clip(s, 260)))
            if not pos:
                silent.append((qname, bucket))
        if len(seen_pos) > 1:
            conflicts.append((qname, sorted(seen_pos)))
            rep.w('    ==> CONFLICT: %s' % ' / '.join(sorted(seen_pos)))
        elif not any_section:
            rep.w('    ==> quantity never mentioned in the swept sections')
        rep.w()

    if conflicts:
        res = 'FAIL'
        detail = '; '.join('%s carries {%s}' % (q, ', '.join(s))
                           for q, s in conflicts)
    elif silent:
        res = 'UNDECIDED'
        detail = ('%d quantity/section pairs mention the quantity with no '
                  'status word, so equality of status cannot be confirmed'
                  % len(silent))
    else:
        res = 'PASS'
        detail = 'every compared quantity carries one status everywhere'
    rep.w('C-01 verdict: %s -- %s' % (res, detail))
    if silent:
        rep.w('  quantity/section pairs with no status word:')
        for q, b in silent:
            rep.w('    - %s in %s' % (q, b))
    rep.w()
    rep.verdict('C-01', 'same status in every section', res, detail)


# ---- C-02 / C-03 ----------------------------------------------------------

def _banned_phrase(ms, rep, cid, title, exact, loose, why):
    rep.w('-' * 78)
    rep.w('[%s]  %s' % (cid, title))
    rep.w('-' * 78)
    rep.w('why banned: ' + why)
    hits_exact, hits_loose = [], []
    for p in ms.paras:
        for sent in sentences(p.text):
            if re.search(exact, sent, re.I):
                hits_exact.append((p.idx, p.bucket or p.section, sent))
            elif loose and re.search(loose, sent, re.I):
                hits_loose.append((p.idx, p.bucket or p.section, sent))
    for idx, sec, s in hits_exact:
        rep.w('  HIT      p%-5d [%s] "%s"' % (idx, sec, clip(s, 320)))
    for idx, sec, s in hits_loose:
        rep.w('  variant  p%-5d [%s] "%s"' % (idx, sec, clip(s, 320)))
    if hits_exact:
        res, det = 'FAIL', '%d occurrence(s)' % len(hits_exact)
    elif hits_loose:
        res, det = 'UNDECIDED', '%d near-variant(s), no exact hit' % len(hits_loose)
    else:
        res, det = 'PASS', 'phrase absent'
        rep.w('  (no occurrence)')
    rep.w('%s verdict: %s -- %s' % (cid, res, det))
    rep.w()
    rep.verdict(cid, title, res, det)
    return hits_exact


def check_c02(ms, rep):
    _banned_phrase(
        ms, rep, 'C-02', 'no occurrence of "the off-diagonal entries converge"',
        exact=r'off-?diagonal entries converge',
        loose=r'off-?diagonal.{0,40}converge',
        why=('the off-diagonal series converges only conditionally and only '
             'for the asymptotic kernel; stated flatly it reads as a property '
             'of the operator'))


def check_c03(ms, rep):
    _banned_phrase(
        ms, rep, 'C-03', 'no occurrence of "more diagonally dominant"',
        exact=r'more diagonally dominant',
        loose=r'diagonal(?:ly)? dominan\w*',
        why=('a matrix whose diagonal grows without bound is not "more '
             'diagonally dominant"; the statement converts a divergence into '
             'a virtue'))


# ---- C-04 -----------------------------------------------------------------

RATIO_SYMBOLIC = re.compile(
    r'\u0394\s*_\{k\s*\+\s*1\}\s*/\s*\u0394\s*_\{k\}\s*(?:\u2192|->)\s*1(?![\d.])')
RATIO_PROSE = re.compile(
    r'(?:ratio of successive increments|increment ratios|successive-increment '
    r'ratios|\u0394\s*_\{k)', re.I)
UNITY_PROSE = re.compile(
    r'(?:\u2192|->|towards|toward|to)\s*(?:unity|1)\b|climb\w* .{0,25}unity|'
    r'rise\w* .{0,25}unity', re.I)
ALPHA_ONE_STRICT = re.compile(
    r'\u03b1\s*=\s*1(?![\d.])|\balpha\s*=\s*1(?![\d.])|'
    r'order\s+exactly\s+(?:one|1)\b|of order one\b(?!\s+or)', re.I)
ALPHA_ONE_WEAK = re.compile(
    r'order one or higher|order\s+(?:one|1)\b|\u03b1\s*\u2265\s*1|'
    r'\u03b1\s*>=\s*1|symbol of order', re.I)


def check_c04(ms, rep):
    cid = 'C-04'
    title = ('"Delta_{k+1}/Delta_k -> 1" only with the alpha = 1 clause')
    rep.w('-' * 78)
    rep.w('[%s]  %s' % (cid, title))
    rep.w('-' * 78)
    rep.w('The ratio tends to 1 for a logarithmic tail, i.e. for a symbol of '
          'order alpha = 1 exactly.')
    rep.w('Accompanying clause must be in the same paragraph and must be '
          'strict (alpha = 1 / "of order one" not followed by "or").')
    rep.w()
    strict_hits, prose_hits = [], []
    for p in ms.paras:
        sents = sentences(p.text)
        has_strict = bool(ALPHA_ONE_STRICT.search(p.text))
        weak = ALPHA_ONE_WEAK.findall(p.text)
        for sent in sents:
            if RATIO_SYMBOLIC.search(sent):
                strict_hits.append((p.idx, p.bucket or p.section, sent,
                                    has_strict, weak))
            elif RATIO_PROSE.search(sent) and UNITY_PROSE.search(sent):
                prose_hits.append((p.idx, p.bucket or p.section, sent,
                                   has_strict, weak))
    unaccompanied = []
    for idx, sec, s, ok, weak in strict_hits:
        rep.w('  SYMBOLIC p%-5d [%s] accompanied=%s' % (idx, sec, ok))
        rep.w('           "%s"' % clip(s, 320))
        if weak:
            rep.w('           nearby weaker clause(s): %s' % sorted(set(weak)))
        if not ok:
            unaccompanied.append((idx, s, weak))
    for idx, sec, s, ok, weak in prose_hits:
        rep.w('  PROSE    p%-5d [%s] accompanied=%s  "%s"'
              % (idx, sec, ok, clip(s, 260)))
        if weak:
            rep.w('           nearby weaker clause(s): %s' % sorted(set(weak)))
    if not strict_hits and not prose_hits:
        rep.w('  (no occurrence)')
        res, det = 'PASS', 'the statement does not appear'
    elif not unaccompanied:
        res, det = 'PASS', 'every symbolic occurrence carries the alpha = 1 clause'
    elif any(w for _, _, w in unaccompanied):
        res = 'UNDECIDED'
        det = ('%d occurrence(s) carry only a weaker clause '
               '("order one or higher" / "alpha >= 1"); an editor must decide '
               'whether that is the alpha = 1 qualification'
               % len(unaccompanied))
    else:
        res = 'FAIL'
        det = '%d occurrence(s) with no order clause at all' % len(unaccompanied)
    rep.w('%s verdict: %s -- %s' % (cid, res, det))
    rep.w()
    rep.verdict(cid, title, res, det)


# ---- C-05 -----------------------------------------------------------------

S6_QUANTITIES = (r'flux linkage|synchronous inductance|self-?inductance|'
                 r'mutual|gap flux density|slot sideband|sideband|cogging|'
                 r'magnet loss|iron loss|torque|efficiency|standstill|'
                 r'machine quantity|back-electromotive')
BOUND_WORDS = r'\bbounds?\b|\bbounded\b|\bbounding\b|\bmajoration\b|\blimits?\b'


def check_c05(ms, rep):
    cid, title = 'C-05', 'equation (32) is not invoked as a bound on a Section 6 quantity'
    rep.w('-' * 78)
    rep.w('[%s]  %s' % (cid, title))
    rep.w('-' * 78)
    offenders, benign = [], []
    for p in ms.paras:
        for sent in sentences(p.text):
            if '(32)' not in sent:
                continue
            is_eq = bool(re.search(r'\(32\)\s*$', sent.strip()))
            has_q = re.search(S6_QUANTITIES, sent, re.I)
            has_b = re.search(BOUND_WORDS, sent, re.I)
            neg = _NEG_RE.search(sent)
            rec = (p.idx, p.section, sent, bool(has_q), bool(has_b), bool(neg))
            if is_eq:
                benign.append(rec + ('the display of (32) itself',))
            elif has_q and has_b and not neg:
                offenders.append(rec + ('quantity + bound, not negated',))
            elif has_q and has_b and neg:
                benign.append(rec + ('quantity + bound, but negated (a denial)',))
            else:
                benign.append(rec + ('no quantity+bound pairing',))
    for rec in offenders:
        rep.w('  OFFENDER p%-5d [S%s] %s' % (rec[0], rec[1], rec[6]))
        rep.w('           "%s"' % clip(rec[2], 340))
    rep.w('  every sentence carrying "(32)", for the record:')
    for rec in benign:
        rep.w('    p%-5d [S%s] (%s)' % (rec[0], rec[1], rec[6]))
        rep.w('           "%s"' % clip(rec[2], 300))
    if offenders:
        res, det = 'FAIL', '%d sentence(s) invoke (32) as a bound' % len(offenders)
    else:
        res, det = 'PASS', ('%d sentence(s) carry (32); none pairs it with a '
                            'Section 6 quantity in a positive bound claim'
                            % (len(benign)))
    rep.w('%s verdict: %s -- %s' % (cid, res, det))
    rep.w()
    rep.verdict(cid, title, res, det)


# ---- C-06 -----------------------------------------------------------------

ATOM = re.compile(r'([A-Za-z\u0370-\u03ff])_\{([^{}]*)\}|'
                  r'([A-Za-z\u0370-\u03ff])(?![A-Za-z_])')

# Symbols a previous manual audit established as used but missing from the
# Nomenclature of the submitted version.  They are a regression guard: the
# check must keep failing while any of them is used and undeclared.
WATCHLIST = ['B_{g1}', 'L_{d}', 'L_{a}', 'A_{ss}', '\u0393', '\u03ba_{n}',
             'b_{0}', '\u03c4_{s}', 'h_{m}', '\u03bc_{I}', 'n_{st}', 'n_{y}']

# tokens that a linearised OMML string yields but that are not symbols
FUNCTION_WORDS = set('''sin cos tan cot sinh cosh tanh coth ln log exp lim
max min sup inf diag rank arctan sqrt det mod'''.split())


def _atoms(expr):
    """Yield normalised symbol keys from one linearised math expression."""
    # strip function names first so that "cosh" does not yield c, o, s, h
    e = expr
    for f in sorted(FUNCTION_WORDS, key=len, reverse=True):
        e = re.sub(r'\b%s\b' % f, ' ', e)
    e = e.replace('\u2211', ' ')
    out = []
    for m in ATOM.finditer(e):
        if m.group(1):
            out.append('%s_{%s}' % (m.group(1), m.group(2)))
        elif m.group(3):
            out.append(m.group(3))
    return out


def check_c06(ms, rep):
    cid, title = 'C-06', 'Nomenclature is complete in both directions'
    rep.w('-' * 78)
    rep.w('[%s]  %s' % (cid, title))
    rep.w('-' * 78)

    # locate the nomenclature table
    nom_first_col, nom_all = [], set()
    nom_para_ids = set()
    found = False
    for rows in ms.tables:
        flat = []
        for cells in rows:
            for c in cells:
                flat.extend(ms.cell_paras(c))
        if any(p.text.strip().lower() == 'nomenclature' for p in flat):
            found = True
            for p in flat:
                nom_para_ids.add(p.idx)
                for e in p.maths:
                    nom_all.update(_atoms(e))
            for cells in rows:
                if not cells:
                    continue
                for p in ms.cell_paras(cells[0]):
                    for e in p.maths:
                        for a in _atoms(e):
                            if a not in nom_first_col:
                                nom_first_col.append(a)
            break
    if not found:
        rep.w('  Nomenclature table not found in the document.')
        rep.verdict(cid, title, 'FAIL', 'no Nomenclature table')
        return

    rep.w('  Nomenclature table found: %d entries in the symbol column, '
          '%d distinct atoms over both columns.'
          % (len(nom_first_col), len(nom_all)))
    rep.w('  declared (symbol column): %s' % ', '.join(nom_first_col))
    rep.w()

    # everything used outside the nomenclature
    used = defaultdict(list)
    for p in ms.paras:
        if p.idx in nom_para_ids:
            continue
        for e in p.maths:
            for a in _atoms(e):
                if len(used[a]) < 8:
                    used[a].append(p.idx)
                elif used[a][-1] != p.idx:
                    pass
    used_keys = set(used)

    # (a) declared but never used
    unused = [s for s in nom_first_col if s not in used_keys]
    rep.w('  (a) declared in the Nomenclature but never used elsewhere:')
    if unused:
        for s in unused:
            rep.w('        MISSING USE  %s' % s)
    else:
        rep.w('        none')
    rep.w()

    # (b) the regression watchlist
    rep.w('  (b) watchlist (symbols a manual audit established as undeclared):')
    wl_missing = []
    for s in WATCHLIST:
        n = len(used.get(s, []))
        declared = s in nom_all
        if n and not declared:
            wl_missing.append(s)
            rep.w('        UNDECLARED   %-10s used in paragraphs %s'
                  % (s, used.get(s, [])[:8]))
        elif not n:
            rep.w('        not used     %-10s (watchlist entry no longer '
                  'appears in the text)' % s)
        else:
            rep.w('        ok           %-10s declared and used' % s)
    rep.w()

    # (c) the general sweep, reported for triage
    rep.w('  (c) other symbols used but not declared (candidates for triage):')
    primary, variants = [], []
    for a in sorted(used_keys):
        if a in nom_all or a in WATCHLIST:
            continue
        base = a.split('_{')[0]
        if '_{' not in a and len(a) == 1 and a.isascii():
            continue          # bare Latin letter: too noisy to decide
        if base in nom_all or base in nom_first_col:
            variants.append(a)
        else:
            primary.append(a)
    rep.w('      primary candidates (base symbol itself undeclared), %d:'
          % len(primary))
    for a in primary:
        rep.w('        %-12s paragraphs %s' % (a, used[a][:6]))
    rep.w('      index variants of a declared base, %d:' % len(variants))
    rep.w('        ' + ', '.join(variants))
    rep.w()

    if unused or wl_missing:
        res = 'FAIL'
        det = ('%d declared-but-unused, %d watchlist symbols used but '
               'undeclared' % (len(unused), len(wl_missing)))
    elif primary:
        res = 'UNDECIDED'
        det = ('%d further symbols are used and not declared; each needs an '
               'editorial decision' % len(primary))
    else:
        res, det = 'PASS', 'nomenclature complete in both directions'
    rep.w('%s verdict: %s -- %s' % (cid, res, det))
    if wl_missing:
        rep.w('  watchlist symbols still undeclared: ' + ', '.join(wl_missing))
    rep.w()
    rep.verdict(cid, title, res, det)


# ---- C-07 -----------------------------------------------------------------

NU8 = re.compile(r'\u03bd\s*=\s*8|first slot sideband|first sideband', re.I)
COMPARATIVE = re.compile(
    r'\bcloser\b|\bclosest\b|\bnearer\b|\bnearest\b|\bbetter\b|\bworse\b|'
    r'\bthe better model\b|\bindistinguishable\b|\bwins\b|\boutperform\w*\b|'
    r'\bimprove\w*\b|\bagainst\b.{0,40}\bfor the operator\b|'
    r'\bthe weaker\b|\bfarther\b|\bfurther off\b', re.I)
DISPERSION_26 = re.compile(
    r'26\.1\d*|\b26\s*%|\u00b1\s*26|26\.126|28\.5|dispersion the tiling', re.I)


def check_c07(ms, rep):
    cid = 'C-07'
    title = 'comparative conclusions on the nu = 8 column carry the +-26 % dispersion'
    rep.w('-' * 78)
    rep.w('[%s]  %s' % (cid, title))
    rep.w('-' * 78)
    rep.w('Table 7(b): the tiling dispersion of the first slot sideband under '
          'the locked chain is 26.126 % (p.c.) and 28.514 % (hat).')
    rep.w('A comparative conclusion on that column without the dispersion is '
          'a conclusion drawn inside the noise.')
    rep.w()
    bad, ok = [], []
    for p in ms.paras:
        para_has_disp = bool(DISPERSION_26.search(p.text))
        for sent in sentences(p.text):
            if not NU8.search(sent) or not COMPARATIVE.search(sent):
                continue
            here = bool(DISPERSION_26.search(sent))
            rec = (p.idx, p.bucket or p.section, sent, here, para_has_disp)
            (ok if (here or para_has_disp) else bad).append(rec)
    for idx, sec, s, here, near in bad:
        rep.w('  UNACCOMPANIED p%-5d [%s]' % (idx, sec))
        rep.w('                "%s"' % clip(s, 340))
    rep.w('  accompanied (dispersion in the same sentence or paragraph):')
    for idx, sec, s, here, near in ok:
        rep.w('    p%-5d [%s] in-sentence=%s  "%s"'
              % (idx, sec, here, clip(s, 240)))
    if not bad and not ok:
        rep.w('  (no comparative conclusion found on this column)')
        res, det = 'UNDECIDED', 'no sentence matched; the detector may be blind'
    elif bad:
        res, det = 'FAIL', '%d comparative conclusion(s) without the dispersion' % len(bad)
    else:
        res, det = 'PASS', 'all %d comparative conclusions carry it' % len(ok)
    rep.w('%s verdict: %s -- %s' % (cid, res, det))
    rep.w()
    rep.verdict(cid, title, res, det)


# ---- C-08 -----------------------------------------------------------------

MESH_RE = re.compile(r'\bmesh(ed)?\b|\brefined mesh\b|\bmesh model\b|'
                     r'\bmesh column\b|\bpolar mesh\b', re.I)
LUMPED_RE = re.compile(r'\blumped\b|\bone-branch-per-tooth\b|\bcoarser model\b',
                       re.I)
ACCURACY_CLAIM = re.compile(
    r'more accurate|more accurately|better\b|superior|outperform\w*|'
    r'improve\w* (?:on|upon)|closer to the reference|is closer\b|'
    r'the closest\b|same ordering|reaches tenths of a percent', re.I)
STRONG_CLAIM = re.compile(r'more accurate|more accurately|superior|'
                          r'outperform\w*', re.I)


def check_c08(ms, rep):
    cid = 'C-08'
    title = 'no sentence claims the mesh iron model is more accurate than the lumped one'
    rep.w('-' * 78)
    rep.w('[%s]  %s' % (cid, title))
    rep.w('-' * 78)
    rep.w('Table 12 is 9 rows to 8 between the two columns, so the ordering is '
          'undecidable and no accuracy claim is available in either direction.')
    rep.w()
    strong, soft = [], []
    for p in ms.paras:
        for sent in sentences(p.text):
            if not (MESH_RE.search(sent) and LUMPED_RE.search(sent)):
                continue
            if not ACCURACY_CLAIM.search(sent):
                continue
            rec = (p.idx, p.bucket or p.section, sent,
                   sorted(set(m.group(0).lower()
                              for m in ACCURACY_CLAIM.finditer(sent))))
            (strong if STRONG_CLAIM.search(sent) else soft).append(rec)
    for idx, sec, s, cues in strong:
        rep.w('  CLAIM     p%-5d [%s] cues=%s' % (idx, sec, cues))
        rep.w('            "%s"' % clip(s, 340))
    for idx, sec, s, cues in soft:
        rep.w('  CANDIDATE p%-5d [%s] cues=%s' % (idx, sec, cues))
        rep.w('            "%s"' % clip(s, 340))
    if strong:
        res, det = 'FAIL', '%d explicit accuracy claim(s)' % len(strong)
    elif soft:
        res = 'UNDECIDED'
        det = ('%d comparative sentence(s) put the two models side by side '
               'with a comparative cue; each needs an editorial reading'
               % len(soft))
    else:
        rep.w('  (no comparative sentence found)')
        res, det = 'PASS', 'no comparative accuracy claim between the two models'
    rep.w('%s verdict: %s -- %s' % (cid, res, det))
    rep.w()
    rep.verdict(cid, title, res, det)


# ---- C-09 -----------------------------------------------------------------

CREDIT_ROLES = ['conceptualis', 'conceptualiz', 'methodology', 'software',
                'validation', 'formal analysis', 'investigation', 'resources',
                'data curation', 'writing', 'visualis', 'visualiz',
                'supervision', 'project administration', 'funding acquisition']


def check_c09(ms, rep, cite_numbers, n_refs, ref_numfmt):
    cid, title = 'C-09', 'front and back matter required by the journal'
    rep.w('-' * 78)
    rep.w('[%s]  %s' % (cid, title))
    rep.w('-' * 78)
    subs = []

    abstract = ' '.join(p.text for p in ms.abstract_paras)
    words = re.findall(r"[A-Za-z0-9\u00c0-\u024f][\w'\u2019\u2013-]*", abstract)
    n_words = len(words)
    subs.append(('abstract <= 250 words', n_words <= 250,
                 '%d words' % n_words))

    kw = [p for p in ms.paras[:40]
          if re.match(r'^\s*key\s?words\b', p.text, re.I)]
    subs.append(('keywords present', bool(kw),
                 clip(kw[0].text, 150) if kw else 'not found'))

    joined = '\n'.join(p.text for p in ms.paras)
    has_contrib = bool(re.search(r'author contributions|CRediT', joined, re.I))
    roles = [r for r in CREDIT_ROLES if r in joined.lower()]
    subs.append(('CRediT statement present', has_contrib and len(roles) >= 4,
                 'heading=%s, %d CRediT role words found' % (has_contrib, len(roles))))

    ci = [p for p in ms.paras
          if re.search(r'competing interest|conflict of interest', p.text, re.I)]
    subs.append(('Declaration of competing interest', bool(ci),
                 clip(ci[-1].text, 160) if ci else 'not found'))

    author_year = re.findall(r'\((?:[A-Z][A-Za-z\u00c0-\u024f\'\u2019-]+'
                             r'(?: et al\.)?,\s*\d{4})\)', joined)
    numbered_list = ref_numfmt in ('decimal', 'ordinal')
    contiguous = (sorted(cite_numbers) == list(range(1, max(cite_numbers) + 1))
                  if cite_numbers else False)
    style_ok = (not author_year) and numbered_list and contiguous and \
               (n_refs == max(cite_numbers) if cite_numbers else False)
    subs.append(('Elsevier numbered reference style', style_ok,
                 'in-text bracketed numerals=%d distinct, author-year hits=%d, '
                 'reference list=%d entries in a %s numbered list, '
                 'cited range contiguous=%s'
                 % (len(cite_numbers), len(author_year), n_refs,
                    ref_numfmt or 'plain', contiguous)))

    for name, ok, det in subs:
        rep.w('  %-38s %-4s  %s' % (name, 'PASS' if ok else 'FAIL', det))
    rep.w('  note: numbering, bracketing and completeness of the reference '
          'list are decidable here; the punctuation of an individual entry is '
          'not, and is not claimed.')
    bad = [n for n, ok, _ in subs if not ok]
    res = 'PASS' if not bad else 'FAIL'
    det = 'all five sub-checks pass' if not bad else 'failing: ' + '; '.join(bad)
    rep.w('%s verdict: %s -- %s' % (cid, res, det))
    rep.w()
    rep.verdict(cid, title, res, det)


# ---- C-10 -----------------------------------------------------------------

EQ_TAG = re.compile(r'\((\d{1,3})\)\s*$')
CITE_GROUP = re.compile(r'\[([0-9][0-9,;\s\u2013\u2014-]*)\]')
SEC_REF = re.compile(r'\bSections?\s+((?:\d+(?:\.\d+)*)'
                     r'(?:\s*(?:,|and|to|\u2013|-)\s*\d+(?:\.\d+)*)*)')
FIG_REF = re.compile(r'\bFigs?\.?\s+((?:\d+(?:\([a-z]\))?)'
                     r'(?:\s*(?:,|and|to|\u2013|-)\s*\d+(?:\([a-z]\))?)*)|'
                     r'\bFigures?\s+((?:\d+)(?:\s*(?:,|and|to|\u2013|-)\s*\d+)*)')
TAB_REF = re.compile(r'\bTables?\s+((?:\d+(?:\([a-z]\))?)'
                     r'(?:\s*(?:,|and|to|\u2013|-)\s*\d+(?:\([a-z]\))?)*)')
EQ_REF = re.compile(r'\((\d{1,3})\)')


def _expand_numbers(s, dash_is_range=True):
    """Expand '1, 2 and 3' / '13-15' into a list of numbers (as strings)."""
    s = s.replace('\u2013', '-').replace('\u2014', '-')
    out = []
    tokens = re.split(r'\s*(?:,|and|to)\s*', s)
    for t in tokens:
        t = t.strip()
        if not t:
            continue
        t = re.sub(r'\([a-z]\)', '', t).strip()
        m = re.fullmatch(r'(\d+(?:\.\d+)*)\s*-\s*(\d+(?:\.\d+)*)', t)
        if m and dash_is_range and '.' not in m.group(1) and '.' not in m.group(2):
            a, b = int(m.group(1)), int(m.group(2))
            if 0 < b - a < 40:
                out.extend(str(x) for x in range(a, b + 1))
                continue
        for piece in re.split(r'\s*-\s*', t):
            piece = piece.strip()
            if re.fullmatch(r'\d+(?:\.\d+)*', piece):
                out.append(piece)
    return out


def collect_numbering(ms):
    """Return the declared and called sequences for equations, tables,
    figures, references and sections."""
    eqs = []           # (number, paragraph)
    for p in ms.paras:
        if not p.maths:
            continue
        tail = p.text.rstrip().rstrip('}').rstrip()
        m = EQ_TAG.search(tail)
        if m:
            eqs.append((int(m.group(1)), p.idx))
    eq_nums = [n for n, _ in eqs]

    tables, figures = [], []
    for p in ms.paras:
        m = ms.CAPTION_RE.match(p.text)
        if m:
            if m.group(1).lower().startswith('tab'):
                tables.append((int(m.group(2)), p.idx))
            else:
                figures.append((int(m.group(2)), p.idx))

    caption_ids = set(i for _, i in tables) | set(i for _, i in figures)
    eq_ids = set(i for _, i in eqs)

    called_eq, called_tab, called_fig = defaultdict(list), defaultdict(list), defaultdict(list)
    cited, sec_refs = defaultdict(list), defaultdict(list)

    max_eq = max(eq_nums) if eq_nums else 0
    for p in ms.paras:
        t = p.text
        if p.idx not in eq_ids:
            for m in EQ_REF.finditer(t):
                n = int(m.group(1))
                if 1 <= n <= max_eq:
                    called_eq[n].append(p.idx)
        for m in TAB_REF.finditer(t):
            for n in _expand_numbers(m.group(1)):
                if p.idx in caption_ids and re.match(r'^\s*Table\s+%s\s*[.:]' % n, t):
                    continue
                called_tab[int(n)].append(p.idx)
        for m in FIG_REF.finditer(t):
            grp = m.group(1) or m.group(2) or ''
            for n in _expand_numbers(grp):
                if p.idx in caption_ids and \
                        re.match(r'^\s*(Figure|Fig\.)\s+%s\s*[.:]' % n, t):
                    continue
                called_fig[int(n)].append(p.idx)
        for m in CITE_GROUP.finditer(t):
            for n in _expand_numbers(m.group(1)):
                if '.' not in n:
                    cited[int(n)].append(p.idx)
        for m in SEC_REF.finditer(t):
            for n in _expand_numbers(m.group(1), dash_is_range=False):
                sec_refs[n].append(p.idx)

    return dict(eqs=eqs, eq_nums=eq_nums, tables=tables, figures=figures,
                called_eq=called_eq, called_tab=called_tab,
                called_fig=called_fig, cited=cited, sec_refs=sec_refs)


def reference_entries(ms):
    """Reference-list paragraphs: the numbered list that follows the
    'References' heading."""
    out, seen_head, numid = [], False, None
    for p in ms.paras:
        if p.style == 'Heading1' and p.text.strip().lower() == 'references':
            seen_head = True
            continue
        if not seen_head:
            continue
        if p.style == 'Heading1':
            break
        if p.text.strip():
            out.append(p)
            if numid is None:
                numid = p.numid
    return out, numid


def numfmt_of(ms, numid):
    if ms.numbering is None or numid is None:
        return None
    abs_id = None
    for num in ms.numbering.findall(W + 'num'):
        if num.get(W + 'numId') == numid:
            a = num.find(W + 'abstractNumId')
            if a is not None:
                abs_id = a.get(W + 'val')
    if abs_id is None:
        return None
    for an in ms.numbering.findall(W + 'abstractNum'):
        if an.get(W + 'abstractNumId') == abs_id:
            lvl = an.find(W + 'lvl')
            if lvl is not None:
                f = lvl.find(W + 'numFmt')
                if f is not None:
                    return f.get(W + 'val')
    return None


def check_c10(ms, rep, num, n_refs, ref_numfmt):
    cid = 'C-10'
    title = 'equations, tables, figures, references and section cross-references'
    rep.w('-' * 78)
    rep.w('[%s]  %s' % (cid, title))
    rep.w('-' * 78)
    problems = []

    def seq_report(label, pairs, called, expect_all_called=True):
        nums = [n for n, _ in pairs]
        arr = np.array(sorted(nums), dtype=int) if nums else np.array([], dtype=int)
        dups = sorted(set(int(x) for x in arr[:-1][arr[1:] == arr[:-1]])) if arr.size > 1 else []
        top = int(arr.max()) if arr.size else 0
        missing = [n for n in range(1, top + 1) if n not in set(nums)]
        never = [n for n in sorted(set(nums)) if not called.get(n)]
        rep.w('  %s: %d declared, numbered 1..%d' % (label, len(nums), top))
        if missing:
            rep.w('      MISSING from the sequence: %s' % missing)
            problems.append('%s missing %s' % (label, missing))
        if dups:
            rep.w('      DUPLICATE numbers: %s' % dups)
            problems.append('%s duplicated %s' % (label, dups))
        if expect_all_called:
            if never:
                rep.w('      NEVER CALLED outside their own caption: %s' % never)
                problems.append('%s never called: %s' % (label, never))
            else:
                rep.w('      every item is called in the text')
        return top

    top_eq = seq_report('equations', num['eqs'], num['called_eq'])
    top_tab = seq_report('tables', num['tables'], num['called_tab'])
    top_fig = seq_report('figures', num['figures'], num['called_fig'])

    cited = num['cited']
    if cited:
        top_ref = max(cited)
        gaps = [n for n in range(1, top_ref + 1) if n not in cited]
        rep.w('  references: %d entries in the list (%s numbering), '
              '%d distinct numbers cited, highest cited [%d]'
              % (n_refs, ref_numfmt or 'plain', len(cited), top_ref))
        if gaps:
            rep.w('      NEVER CITED: %s' % gaps)
            problems.append('references never cited: %s' % gaps)
        if n_refs != top_ref:
            rep.w('      LIST LENGTH %d != highest cited %d' % (n_refs, top_ref))
            problems.append('reference list length %d vs highest cited %d'
                            % (n_refs, top_ref))
        if not gaps and n_refs == top_ref:
            rep.w('      every entry is cited and the list length matches')
    else:
        rep.w('  references: no bracketed citation found')
        problems.append('no citations parsed')

    existing = set(ms.sections)
    dangling = OrderedDict()
    for ref, where in sorted(num['sec_refs'].items()):
        if ref not in existing:
            dangling[ref] = where
    rep.w('  section cross-references: %d distinct targets, %d exist'
          % (len(num['sec_refs']), len(num['sec_refs']) - len(dangling)))
    if dangling:
        for ref, where in dangling.items():
            rep.w('      DANGLING "Section %s" cited in paragraphs %s'
                  % (ref, sorted(set(where))[:8]))
        problems.append('dangling section references: %s' % list(dangling))
    else:
        rep.w('      every "Section x.y" resolves to a heading')
    rep.w('  section map:')
    for k, v in ms.sections.items():
        rep.w('      %-8s %s' % (k, v))

    res = 'PASS' if not problems else 'FAIL'
    det = 'complete sequences, all called, no dangling reference' if not problems \
        else '; '.join(problems)
    rep.w('%s verdict: %s -- %s' % (cid, res, det))
    rep.w()
    rep.verdict(cid, title, res, det)


# --------------------------------------------------------------------------
# 5.  main
# --------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('docx', nargs='?', default=DEFAULT_DOCX,
                    help='manuscript .docx (default: the submitted version)')
    ap.add_argument('--out', default=DEFAULT_OUT,
                    help='transcript file to write')
    args = ap.parse_args(argv)

    if not os.path.isfile(args.docx):
        sys.stderr.write('manuscript not found: %s\n' % args.docx)
        return 2

    ms = Manuscript(args.docx)
    rep = Report()

    rep.w('=' * 78)
    rep.w('B12  CONSISTENCY CHECK  --  terminology and status')
    rep.w('=' * 78)
    rep.w('Date (ISO)   : %s' % ISO_DATE)
    rep.w('Manuscript   : %s' % args.docx)
    rep.w('Size (bytes) : %d' % os.path.getsize(args.docx))
    rep.w('Script       : %s' % os.path.abspath(__file__))
    rep.w('Paragraphs   : %d (%d inside w:tbl)'
          % (len(ms.paras), sum(1 for p in ms.paras if p.in_table)))
    rep.w('OMML expr.   : %d' % sum(len(p.maths) for p in ms.paras))
    rep.w('Sections     : %d numbered headings' % len(ms.sections))
    rep.w('Rule         : a check that cannot be settled mechanically reports '
          'UNDECIDED, never PASS.')
    rep.w()

    num = collect_numbering(ms)
    refs, ref_numid = reference_entries(ms)
    ref_numfmt = numfmt_of(ms, ref_numid)

    check_c01(ms, rep)
    check_c02(ms, rep)
    check_c03(ms, rep)
    check_c04(ms, rep)
    check_c05(ms, rep)
    check_c06(ms, rep)
    check_c07(ms, rep)
    check_c08(ms, rep)
    check_c09(ms, rep, set(num['cited']), len(refs), ref_numfmt)
    check_c10(ms, rep, num, len(refs), ref_numfmt)

    rep.w('=' * 78)
    rep.w('SUMMARY')
    rep.w('=' * 78)
    for cid, (title, res, det) in rep.verdicts.items():
        rep.w('%-6s %-10s %s' % (cid, res, title))
        rep.w('       %s' % det)
    n_fail = sum(1 for _, r, _ in rep.verdicts.values() if r == 'FAIL')
    n_und = sum(1 for _, r, _ in rep.verdicts.values() if r == 'UNDECIDED')
    n_pass = sum(1 for _, r, _ in rep.verdicts.values() if r == 'PASS')
    rep.w()
    rep.w('PASS %d   FAIL %d   UNDECIDED %d   (of %d checks)'
          % (n_pass, n_fail, n_und, len(rep.verdicts)))
    ok = (n_fail == 0 and n_und == 0)
    rep.w('exit code: %d' % (0 if ok else 1))

    out = rep.text()
    d = os.path.dirname(os.path.abspath(args.out))
    if d and not os.path.isdir(d):
        os.makedirs(d)
    with open(args.out, 'w', encoding='utf-8') as f:
        f.write(out)
    sys.stdout.write(out)
    sys.stdout.write('\ntranscript written to %s\n' % os.path.abspath(args.out))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
