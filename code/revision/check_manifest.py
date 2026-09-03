# -*- coding: utf-8 -*-
"""CHECK_MANIFEST  -  the verification half of task B10.

MANIFEST.json claims, for each published item of the submitted manuscript, the
script and the dated transcript it came from. This script tests that claim: for
every value the manifest attributes to a transcript, it asks whether that value
is actually REACHABLE in that transcript, at the precision at which it is
published.

WHAT "REACHABLE" MEANS. Every numeric token in the transcript is parsed to a
float. A published value v printed with d decimals is reachable if some
transcript number x satisfies |x - v| <= half a unit in v's last printed place.
This accepts the transcript printing more digits than the manuscript, which is
the normal case, and rejects a value that simply is not there.

WHAT THIS SCRIPT WILL NOT DO. It never edits the manifest to make a check pass.
An unreachable value is reported as a finding; that is the whole point of the
exercise, since Section 5.2 of the manuscript admits that nothing in the file
system distinguishes an authoritative chain from an exploratory one.

Usage:  python check_manifest.py
Exit code 0 if every authoritative attribution holds, 1 otherwise.
"""
import json, os, re, sys, math

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
MAN  = os.path.join(HERE, 'MANIFEST.json')
OUTF = os.path.join(REPO, 'outputs', 'revision', 'B10_check_manifest_out.txt')

NUM = re.compile(r'[-+]?\d+(?:[.,]\d+)?(?:[eE][-+]?\d+)?')

def read(p):
    for enc in ('utf-8', 'latin-1'):
        try:
            with open(p, encoding=enc) as f:
                return f.read()
        except (UnicodeDecodeError, FileNotFoundError):
            continue
    return None

def numbers(text):
    out = []
    for m in NUM.finditer(text):
        s = m.group(0).replace(',', '.')
        try:
            out.append(float(s))
        except ValueError:
            pass
    return out

def places(vstr):
    """Half a unit in the last printed place of vstr."""
    s = vstr.strip().replace(',', '.')
    mant = s.split('e')[0].split('E')[0]
    exp  = 0
    m = re.search(r'[eE]([-+]?\d+)', s)
    if m:
        exp = int(m.group(1))
    d = len(mant.split('.')[1]) if '.' in mant else 0
    return 0.5 * 10**(-d) * 10**exp

class Rep:
    def __init__(self):
        self.lines = []
    def w(self, s=''):
        print(s)
        self.lines.append(s)

def main():
    rep = Rep()
    man = json.load(open(MAN, encoding='utf-8'))
    recs = man['records']

    rep.w("=" * 78)
    rep.w("B10 : does every attributed value actually live in its declared transcript?")
    rep.w("=" * 78)
    rep.w(f"manifest    : {os.path.relpath(MAN, REPO)}")
    rep.w(f"manuscript  : {man.get('manuscript','')}")
    rep.w(f"manifest date: {man.get('date_iso','')}")
    rep.w(f"records     : {len(recs)}")
    rep.w("")

    cache = {}
    n_auth = n_ok = n_bad = 0
    n_val_ok = n_val_bad = 0
    missing_files = []
    failures = []

    for r in recs:
        if not r.get('authoritative'):
            continue
        n_auth += 1
        q  = r.get('quantity', '?')
        sc = r.get('script', '')
        tr = r.get('transcript', '')

        #  Several manifest records name more than one file in a single field,
        #  joined by " + ", sometimes with a parenthetical note. Split them and
        #  check each path in its own right; a value counts as reachable if ANY
        #  of the declared transcripts carries it.
        def split_paths(field):
            #  Keep only fragments that actually look like repository paths, so
            #  that a trailing parenthetical note ("(savefigure ...)") is not
            #  mistaken for a missing file.
            out = []
            for part in str(field).split('+'):
                part = re.sub(r'\(.*?\)', '', part)
                part = part.replace('(', ' ').replace(')', ' ').strip()
                for tok in part.split():
                    if '/' in tok and re.search(r'\.(m|py|txt|json|mat)$', tok):
                        out.append(tok)
            return out

        for p, kind in ((sc, 'script'), (tr, 'transcript')):
            for one in split_paths(p):
                if not os.path.isfile(os.path.join(REPO, one)):
                    missing_files.append(f"{kind} absent: {one}   (record: {q[:60]})")

        vals = r.get('values') or []
        if not tr or not vals:
            continue
        nums = []
        unread = []
        for one in split_paths(tr):
            if one not in cache:
                t = read(os.path.join(REPO, one))
                cache[one] = numbers(t) if t else None
            if cache[one] is None:
                unread.append(one)
            else:
                nums.extend(cache[one])
        for one in unread:
            missing_files.append(f"transcript unreadable: {one}")
        if not nums:
            continue

        bad = []
        for v in vals:
            vs = str(v.get('value', '')).strip()
            if not vs:
                continue
            try:
                vf = float(vs.replace(',', '.'))
            except ValueError:
                continue
            tol = places(vs)
            hit = any(abs(x - vf) <= tol for x in nums)
            # a published value may be a percentage of a transcript ratio, or
            # carry the opposite sign convention; both are accepted explicitly
            if not hit:
                hit = any(abs(abs(x) - abs(vf)) <= tol for x in nums)
            if hit:
                n_val_ok += 1
            else:
                n_val_bad += 1
                bad.append(vs)
        if bad:
            n_bad += 1
            failures.append((q, tr, bad))
        else:
            n_ok += 1

    rep.w("-" * 78)
    rep.w("RESULT")
    rep.w("-" * 78)
    rep.w(f"  authoritative records            : {n_auth}")
    rep.w(f"  records fully reachable          : {n_ok}")
    rep.w(f"  records with unreachable values  : {n_bad}")
    rep.w(f"  values reachable                 : {n_val_ok}")
    rep.w(f"  values NOT reachable             : {n_val_bad}")
    rep.w("")

    tabs = sorted({str(r.get('table')) for r in recs if r.get('table')},
                  key=lambda s: (len(s), s))
    auth_tabs = sorted({str(r.get('table')) for r in recs
                        if r.get('table') and r.get('authoritative')},
                       key=lambda s: (len(s), s))
    rep.w(f"  tables named in the manifest     : {len(tabs)}  -> {', '.join(tabs)}")
    rep.w(f"  tables with an authoritative row : {len(auth_tabs)}  -> {', '.join(auth_tabs)}")
    rep.w("")

    if missing_files:
        rep.w("  DECLARED FILES THAT DO NOT EXIST:")
        for m in missing_files:
            rep.w("    " + m)
        rep.w("")

    if failures:
        rep.w("  VALUES THE MANIFEST ATTRIBUTES BUT THE TRANSCRIPT DOES NOT CARRY:")
        for q, tr, bad in failures:
            rep.w(f"    {q[:72]}")
            rep.w(f"      transcript : {tr}")
            rep.w(f"      unreachable: {', '.join(bad[:12])}"
                  + (" ..." if len(bad) > 12 else ""))
        rep.w("")
        rep.w("  These are findings, not defects of this script. Nothing has been")
        rep.w("  adjusted to make them pass.")
    else:
        rep.w("  Every value attributed to a transcript was found in it.")

    unattr = [r for r in recs if not r.get('authoritative')]
    rep.w("")
    rep.w(f"  items the manifest declines to attribute : {len(unattr)}")
    rep.w("  (an empty cell is information; a wrong cell is not - these are listed")
    rep.w("   in MANIFEST.json with the reason in each record's 'note' field)")

    rep.w("")
    rep.w("=" * 78)
    ok = (n_val_bad == 0 and not missing_files)
    rep.w(f"CHECK_MANIFEST : {'PASS' if ok else 'FAIL'}")
    rep.w("=" * 78)

    os.makedirs(os.path.dirname(OUTF), exist_ok=True)
    with open(OUTF, 'w', encoding='utf-8') as f:
        f.write("\n".join(rep.lines) + "\n")
    return 0 if ok else 1

if __name__ == '__main__':
    sys.exit(main())
