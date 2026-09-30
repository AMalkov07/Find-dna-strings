#!/usr/bin/env python3
"""Collect every circle candidate population mode flagged (families >= 90 bp, the ones marked
'<- circle' in the CIRCLE SUMMARY) from each sample's _output_pattern.txt.
Usage: build_units.py W SAMPLES_FILE OUTDIR   -> OUTDIR/units.tsv, OUTDIR/units_info.tsv"""
import re
import sys

W, SF, OUT = sys.argv[1:4]
SAMPLES = [l.strip() for l in open(SF) if l.strip()]


def minper(s):
    for p in range(1, 13):
        if sum(s[i] != s[i + p] for i in range(len(s) - p)) <= 0.03 * len(s):
            return p


rows = []
for s in SAMPLES:
    t = open(f"{W}/{s}/{s}_output_pattern.txt").read()
    short = re.sub(r"^\d+_", "", s).replace("_with_selection", "").replace("survivor_", "").replace("_repeat", "rep")
    if "REPRESENTATIVE SEQUENCES" not in t:
        continue
    flagged = re.findall(r"^\s+(\d+)\s+\d+/\d+.*<- circle", t.split("CIRCLE SUMMARY", 1)[1], re.M)
    fams = re.findall(r"Family #(\d+)\s+\((\d+) bp representative.*?fam_reads=(\d+)/(\d+).*?\n\s+([ACGT]+)",
                      t.split("REPRESENTATIVE SEQUENCES", 1)[1], re.S)
    got = [(f, l, fr, q) for f, l, fr, _, q in fams if int(l) >= 90]
    assert len(got) == len(flagged), (s, len(got), flagged)
    for f, l, fr, q in got:
        p = minper(q)
        rows.append((f"{short}_f{f}_{l}", s, f"family{f}", q,
                     f"simple repeat ({q[:p]})" if p else "TG telomeric", fr))
with open(f"{OUT}/units.tsv", "w") as f:
    for r in rows:
        f.write("\t".join(r[:4]) + "\n")
with open(f"{OUT}/units_info.tsv", "w") as f:
    f.write("unit\tsample\tfamily\tkind\tfam_reads\n")
    for r in rows:
        f.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[4]}\t{r[5]}\n")
print(f"{len(rows)} candidates from {len(SAMPLES)} samples")
