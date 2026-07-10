#!/usr/bin/env python3
"""
Regenerate population *_output_pattern.txt (and its companion *_pattern_reads.txt)
from the existing *_output_telomeres.fasta, using the CURRENT finder code
(new scoring, length-aware consolidation at vt=0.93, CIRCLE SUMMARY verdict).

Telomere extraction is UNCHANGED by recent edits, so the extracted telomeres are
already current — we skip extraction entirely and re-analyze them directly. No raw
reads needed, no re-extraction, much lighter than a full pipeline rerun.

Usage:  python refresh_from_telomeres.py <population_results_dir>
Reads NSLOTS (SGE) for the worker count if set.
"""
import glob
import os
import re
import sys

from Bio import SeqIO
from analysis.pattern_finder_lcp import analyze_population_reads


def raw_read_count(pattern_path: str, fallback: int) -> int:
    """Carry over the 'Total raw reads' header from the existing pattern file
    (we don't have the raw reads here). Falls back to the telomere count."""
    try:
        with open(pattern_path) as f:
            m = re.search(r"Total raw reads\s*:\s*(\d+)", f.read())
            if m:
                return int(m.group(1))
    except OSError:
        pass
    return fallback


def main(argv) -> int:
    results_dir = argv[0] if argv else "."
    workers = int(os.environ.get("NSLOTS", "0")) or None
    tels = sorted(glob.glob(os.path.join(results_dir, "*_output_telomeres.fasta")))
    if not tels:
        print(f"No *_output_telomeres.fasta found in {results_dir}", file=sys.stderr)
        return 1
    print(f"Refreshing {len(tels)} sample(s) from telomeres  (workers={workers or 'auto'})")
    for tel in tels:
        pat = tel.replace("_output_telomeres.fasta", "_output_pattern.txt")
        seqs = {r.id: str(r.seq) for r in SeqIO.parse(tel, "fasta")}
        if not seqs:
            print(f"  SKIP {os.path.basename(tel)}: no reads", flush=True)
            continue
        traw = raw_read_count(pat, len(seqs))
        print(f"  {os.path.basename(pat)}: {len(seqs)} telomere reads (raw={traw}) ...",
              flush=True)
        with open(pat, "w") as pf:
            analyze_population_reads(seqs, pf, 50, 300,
                                     total_raw_reads=traw, max_workers=workers)
        print("    done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
