#!/usr/bin/env python3
"""
Second-pass read classifier for consolidated circle families.

For each family in a population *_output_pattern.txt (its representative sequence,
from the "CONSOLIDATED CIRCLE FAMILIES — REPRESENTATIVE SEQUENCES" section), scan
every extracted telomere read (<base>_output_telomeres.fasta) and classify it by
how the family's circle appears in it:

  CIRCULARIZED    - the read contains a TANDEM of the circle (>= 2 adjacent copies),
                    i.e. the doubled representative aligns into the read at >= identity.
                    Strongest evidence (rolling-circle signature).
  MULTI-COPY      - the read contains >= 2 copies of the circle that are NOT tandem
                    (dispersed across the read). Medium evidence: stronger than one
                    copy, weaker than a clean tandem.
  SINGLE INSTANCE - the read contains exactly ONE clean copy. Weakest; the first
                    (suffix-array) pass misses these because it needs >=2 in-read repeats.

Matching is rotation-aware (via the doubled representative) and error-tolerant
(edit distance), so nanopore errors don't break it. It is deliberately strict
(full-length, high identity) to avoid pulling in GT-rich look-alikes.

Requires: edlib  (pip install edlib), biopython.

Usage:
    python classify_family_reads.py FILE_output_pattern.txt [more ...]
    python classify_family_reads.py --identity 0.90 --top-families 6 FILE...
Writes <base>_family_reads.txt next to each input.
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Tuple

try:
    import edlib
except ImportError:
    sys.exit("ERROR: this script needs the 'edlib' package (pip install edlib)")
try:
    from Bio import SeqIO
except ImportError:
    sys.exit("ERROR: this script needs biopython (pip install biopython)")

DEFAULT_IDENTITY = 0.98
COPY_CAP = 8  # cap the reported tandem-copy estimate

_FAMILY_RE = re.compile(
    r"Family #(\d+)\s+\((\d+) bp representative.*?\n\s+([ACGTNacgtn]+)",
    re.MULTILINE,
)


def parse_families(text: str) -> List[Tuple[int, str]]:
    """Return [(family_number, representative_seq), ...]."""
    if "REPRESENTATIVE SEQUENCES" in text:
        section = text.split("REPRESENTATIVE SEQUENCES", 1)[1]
    else:
        section = text
    return [(int(n), seq.upper()) for n, _, seq in _FAMILY_RE.findall(section)]


def hw_identity(query: str, target: str) -> float:
    """Best identity placing the FULL query inside target (free target ends)."""
    if not query or len(query) > len(target):
        return 0.0
    d = edlib.align(query, target, mode="HW", task="distance")["editDistance"]
    return 1.0 - d / len(query)


def copy_estimate(rep: str, read: str, identity: float) -> int:
    """Largest k (2..CAP) for which rep*k still aligns into read at >= identity."""
    best = 1
    for k in range(2, COPY_CAP + 1):
        q = rep * k
        if len(q) > len(read):
            break
        if hw_identity(q, read) >= identity:
            best = k
        else:
            break
    return best


def count_occurrences(rep: str, read: str, identity: float, cap: int = 20) -> int:
    """Number of non-overlapping full-length matches of rep in read at >= identity.
    Greedy: repeatedly take the best match, cut it out, re-scan."""
    n = 0
    s = read
    L = len(rep)
    while len(s) >= L and n < cap:
        r = edlib.align(rep, s, mode="HW", task="locations")
        if not r["locations"] or (1 - r["editDistance"] / L) < identity:
            break
        st, en = r["locations"][0]
        s = s[:st] + s[en + 1:]   # remove matched region, look for the next copy
        n += 1
    return n


def classify_read(rep: str, rep2: str, read: str, identity: float) -> Tuple[str, float, int]:
    """Return (category, best_identity, copies).
    category in {circularized, multi, single, ''}:
      circularized - >=2 TANDEM copies (doubled rep aligns; rolling-circle signature)
      multi        - >=2 copies but NOT tandem (dispersed across the read)
      single       - exactly one clean copy
    """
    if len(read) < len(rep):
        return "", 0.0, 0
    # tandem (>=2 adjacent copies) — strongest evidence
    id2 = hw_identity(rep2, read)
    if id2 >= identity:
        return "circularized", id2, copy_estimate(rep, read, identity)
    id1 = hw_identity(rep, read)
    if id1 < identity:
        return "", max(id1, 0.0), 0
    # at least one copy, but no tandem — count how many (dispersed)
    occ = count_occurrences(rep, read, identity)
    if occ >= 2:
        return "multi", id1, occ
    return "single", id1, 1


def process(pattern_path: str, identity: float, top_families: int,
            min_rep_length: int) -> None:
    with open(pattern_path) as f:
        text = f.read()
    families = parse_families(text)
    if not families:
        print(f"  SKIP {pattern_path}: no family representatives found", file=sys.stderr)
        return
    # Drop short families (generic telomere GT-repeat backbone, not distinct
    # circles) — the main lever against listing too many reads.
    families = [(n, s) for (n, s) in families if len(s) >= min_rep_length]
    if top_families:
        families = families[:top_families]
    if not families:
        print(f"  SKIP {pattern_path}: no families >= {min_rep_length} bp", file=sys.stderr)
        return

    tel_path = pattern_path.replace("_output_pattern.txt", "_output_telomeres.fasta")
    if not os.path.isfile(tel_path):
        print(f"  SKIP {pattern_path}: telomeres file not found ({tel_path})", file=sys.stderr)
        return
    reads: Dict[str, str] = {r.id: str(r.seq).upper() for r in SeqIO.parse(tel_path, "fasta")}

    out_path = pattern_path.replace("_output_pattern.txt", "_output_family_reads.txt")
    with open(out_path, "w") as out:
        out.write("READS PER CIRCLE FAMILY (second pass)\n")
        out.write("=" * 90 + "\n")
        out.write(f"  telomere reads scanned : {len(reads)}\n")
        out.write(f"  identity threshold      : {identity:.0%} (full-length, rotation/indel tolerant)\n")
        out.write("  CIRCULARIZED = read has >=2 TANDEM copies (rolling-circle signature; strongest)\n")
        out.write("  MULTI-COPY   = read has >=2 copies but NOT tandem (dispersed; medium)\n")
        out.write("  SINGLE       = read has exactly one clean copy (weakest)\n")
        out.write("=" * 90 + "\n\n")

        for fam_no, rep in families:
            rep2 = rep + rep
            circ: List[Tuple[str, float, int]] = []
            multi: List[Tuple[str, float, int]] = []
            single: List[Tuple[str, float, int]] = []
            for rid, seq in reads.items():
                cat, ident, copies = classify_read(rep, rep2, seq, identity)
                if cat == "circularized":
                    circ.append((rid, ident, copies))
                elif cat == "multi":
                    multi.append((rid, ident, copies))
                elif cat == "single":
                    single.append((rid, ident, copies))
            circ.sort(key=lambda x: (-x[2], -x[1]))
            multi.sort(key=lambda x: (-x[2], -x[1]))
            single.sort(key=lambda x: -x[1])

            out.write("#" * 90 + "\n")
            out.write(f"FAMILY {fam_no}  |  representative {len(rep)} bp  |  "
                      f"circularized: {len(circ)}  |  multi-copy: {len(multi)}  |  single: {len(single)}\n")
            out.write(f"representative: {rep}\n")
            out.write("#" * 90 + "\n")

            out.write(f"\n--- CIRCULARIZED (>=2 tandem copies, >= {identity:.0%} id) : {len(circ)} reads ---\n")
            for rid, ident, copies in circ:
                cp = f"{copies}{'+' if copies >= COPY_CAP else ''}"
                out.write(f">{rid}  copies={cp}  id={ident:.2f}\n{reads[rid]}\n")

            out.write(f"\n--- MULTI-COPY (>=2 non-tandem copies, >= {identity:.0%} id) : {len(multi)} reads ---\n")
            for rid, ident, copies in multi:
                cp = f"{copies}{'+' if copies >= 20 else ''}"
                out.write(f">{rid}  copies={cp}  id={ident:.2f}\n{reads[rid]}\n")

            out.write(f"\n--- SINGLE INSTANCE (1 copy, >= {identity:.0%} id) : {len(single)} reads ---\n")
            for rid, ident, _ in single:
                out.write(f">{rid}  copies=1  id={ident:.2f}\n{reads[rid]}\n")
            out.write("\n")

            print(f"  {os.path.basename(pattern_path)}  Family {fam_no} ({len(rep)}bp): "
                  f"{len(circ)} circularized, {len(multi)} multi-copy, {len(single)} single")
    print(f"  -> {out_path}")


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="population *_output_pattern.txt file(s)")
    ap.add_argument("-i", "--identity", type=float, default=DEFAULT_IDENTITY,
                    help=f"minimum full-length identity to list a read (default {DEFAULT_IDENTITY}). "
                         "Telomeric reads share a ~0.88-0.92 GT-rich baseline, so 0.90 barely beats "
                         "chance; a real circle copy is near read accuracy, so 0.98 is used to require "
                         "a genuine match to the specific circle rather than generic telomere similarity")
    ap.add_argument("-n", "--top-families", type=int, default=0,
                    help="only process the first N families (0 = all)")
    ap.add_argument("-m", "--min-rep-length", type=int, default=100,
                    help="skip families whose representative is shorter than this many bp "
                         "(default 100: drops the generic short telomere-repeat families so "
                         "you don't list most of the reads). Set 0 to process all families.")
    args = ap.parse_args(argv)
    for path in args.files:
        process(path, args.identity, args.top_families, args.min_rep_length)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
