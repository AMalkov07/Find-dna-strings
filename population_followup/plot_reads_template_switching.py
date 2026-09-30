#!/usr/bin/env python3
"""Template-switching graph of every telomeric read carrying one circle unit, one row per read.

Uses Find-dna-strings' own TemplateSwitchingStrategy (analysis, program defaults) and
GraphTemplateSwitching.graph_alignment_analysis (drawing), i.e. the reference-mode template-switching
graph, but each read gets its own right-pointing row (odd chr_index -> sign +1,
y = 17 - (chr_index // 2) / 2), so there is no 32-chromosome-end limit. Reads are drawn in one block
per sample, sorted by copy number.

Each SOURCE is LABEL=FILE, where FILE is either
  - a population-search matches FASTA (<sample>_<unit>_perfect_search_matches.fasta; exact copies), or
  - a population-mode family read list (<sample>_output_family_reads.txt) followed by ':<family number>'
    to take that family's reads (>= 98% identity; use this when exact matching finds few reads).
Headers must carry 'copies=N' (both formats do).

Usage (Find-dna-strings root on PYTHONPATH):
    python population_followup/plot_reads_template_switching.py --unit-seq SEQ | --units-tsv units.tsv --unit NAME
        --out OUT_PREFIX [--title TEXT] LABEL=FILE[:FAMILY] [LABEL=FILE[:FAMILY] ...]
Writes OUT_PREFIX.png and OUT_PREFIX_rows.tsv (row -> read id, copies, length, segments)."""
import argparse
import re
import sys

import matplotlib
matplotlib.use("Agg")

from utils.data_structures import Config, TelomereSequence
from analysis.template_switching_strategy_2 import TemplateSwitchingStrategy
from data_io.graph_template_switching import GraphTemplateSwitching

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("sources", nargs="+", help="LABEL=FILE or LABEL=FAMILY_READS_FILE:FAMILY")
ap.add_argument("--unit-seq", help="circle unit sequence (TG orientation, as population mode reports it)")
ap.add_argument("--units-tsv", help="units.tsv from build_units.py (name, sample, family, sequence)")
ap.add_argument("--unit", help="unit name to take from --units-tsv")
ap.add_argument("--out", required=True, help="output prefix")
ap.add_argument("--title", default="", help="text placed before the standard title")
args = ap.parse_args()

if args.unit_seq:
    PATTERN, UNIT = args.unit_seq.strip().upper(), f"{len(args.unit_seq.strip())} bp unit"
elif args.units_tsv and args.unit:
    PATTERN = next(l.rstrip("\n").split("\t")[3] for l in open(args.units_tsv) if l.split("\t")[0] == args.unit)
    UNIT = args.unit
else:
    sys.exit("give --unit-seq, or --units-tsv with --unit")


def read_source(spec):
    """Return [(read_id, copies, sequence)] from a matches FASTA or a family-reads block."""
    fam = None
    m = re.match(r"(.+_family_reads\.txt):(\d+)$", spec)
    if m:
        spec, fam = m.group(1), m.group(2)
    text = open(spec).read()
    if fam is not None:
        blocks = re.split(r"#{10,}\nFAMILY ", text)
        block = next((b for b in blocks if b.split()[0] == fam), None)
        if block is None:
            sys.exit(f"family {fam} not found in {spec}")
        text = block
    return [(rid, int(c), s) for rid, c, s in
            re.findall(r"^>(\S+)[^\n]*?copies=(\d+)\+?[^\n]*\n([ACGTNacgtn]+)", text, re.M)]


cfg = Config(fasta_file_path="", output_file="", analysis_strategy="template_switching",
             max_ends=0, pattern=PATTERN, maximum_alignment_mutations=12, skip_seeding=True,
             compare_file_path=None, min_pattern_length=50, max_pattern_length=300,
             graph_output=f"{args.out}.png", mutation_lookahead=40, max_deletion_size=3,
             max_insertion_size=3, no_mutations=False)

sections = []
for src in args.sources:
    label, spec = src.split("=", 1)
    recs = sorted(read_source(spec), key=lambda r: (-r[1], -len(r[2])))
    if recs:
        sections.append((label, recs))
if not sections:
    sys.exit(f"{UNIT}: no reads in any source - no graph")

telomeres = [TelomereSequence(survivor_name=label, sequence=s.upper(), chromosome_end_id=rid, analysis=None)
             for label, recs in sections for rid, _, s in recs]
TemplateSwitchingStrategy(telomeres, PATTERN, cfg).execute()

n_rows = sum(len(r) for _, r in sections) + 3 * len(sections)
grapher = GraphTemplateSwitching([], PATTERN, cfg)
grapher.fig.set_size_inches(16, max(4, 0.14 * n_rows))
grapher.ax.clear()
grapher.ax.get_yaxis().set_visible(False)
grapher.ax.get_xaxis().set_visible(False)

LABEL_X = 260
k, it = 0, iter(telomeres)
with open(f"{args.out}_rows.tsv", "w") as tsv:
    tsv.write("row\tsample\tread_id\tcopies\tlength_bp\tn_segments\n")
    for label, recs in sections:
        grapher.ax.text(LABEL_X, 17 - k / 2, f"{label}:  {len(recs)} reads with {UNIT}",
                        ha="right", va="center", fontsize=6, fontweight="bold")
        k += 2
        for row, (rid, copies, s) in enumerate(recs, 1):
            tel = next(it)
            grapher.graph_alignment_analysis(tel, 2 * k + 1)
            grapher.ax.text(LABEL_X, 17 - k / 2, f"{label} #{row}  {rid[:8]}  {copies}x  {len(s)} bp",
                            ha="right", va="center", fontsize=3.2, family="monospace")
            tsv.write(f"{label}#{row}\t{label}\t{rid}\t{copies}\t{len(s)}\t"
                      f"{len(tel.analysis.template_switch_events)}\n")
            k += 1
        k += 1

grapher.ax.set_ylim(17 - k / 2 - 0.5, 17 + 0.75)
grapher.ax.set_xlim(-1100, grapher.max_x + grapher.perfect_arrow_distance * 1.5)
grapher.ax.set_title(f"{args.title + ': ' if args.title else ''}every read with {UNIT} ({len(PATTERN)} bp) - template "
                     "switching (teal = continuous copy of the unit, blue tick = one unit length, "
                     "+/-N = jump within the unit, gold/purple/red = ins/del/mismatch)", fontsize=6)
grapher.fig.savefig(f"{args.out}.png", dpi=300, bbox_inches="tight")
print(f"{UNIT}: {len(telomeres)} reads in {len(sections)} samples -> {args.out}.png")
