#!/usr/bin/env python3
"""Extract telomere ends ONCE (program's population-mode reader, no length gate), then run the
program's own population search for every unit in OUTDIR/units.tsv. Doing all units in one pass
avoids re-reading a multi-GB FASTA once per unit.

Default is exact matching only ('perfect', what the heatmap and timecourse use); add 'alignment'
to MODES to also write the error-tolerant search (it only analyses reads that have >= 1 exact copy).
Output names match main.py -pm -psm: OUTDIR/<sample>_<unit>_<mode>_search.txt (+ _matches.fasta).

Usage (Find-dna-strings root on PYTHONPATH):
    python population_followup/search_units.py FASTA OUTDIR SAMPLE WORKERS [MODES]
    MODES = comma list, default "perfect", e.g. "perfect,alignment"."""
import sys
from utils.data_structures import Config
from data_io.fasta_reader import FastaReader
from analysis import population_search as ps

fa, outdir, sample, w = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
MODES = tuple(sys.argv[5].split(",")) if len(sys.argv) > 5 else ("perfect",)
units = [l.rstrip("\n").split("\t") for l in open(f"{outdir}/units.tsv")]
tel, raw = FastaReader(fa, max_ends=0).parse_fasta_population(threshold=0.75, max_workers=w)
print(f"{sample}: {raw} reads -> {len(tel)} telomere ends (no length gate)")
for name, _, _, seq in units:
    for mode in MODES:
        out = f"{outdir}/{sample}_{name}_{mode}.txt"
        cfg = Config(fasta_file_path=fa, output_file=out, analysis_strategy="template_switching",
                     max_ends=0, pattern=seq, maximum_alignment_mutations=12, skip_seeding=True,
                     compare_file_path=None, min_pattern_length=50, max_pattern_length=300,
                     graph_output=None, mutation_lookahead=40, max_deletion_size=3,
                     max_insertion_size=3, no_mutations=False, population_mode=True,
                     workers=w, population_search_mode=mode, min_telomere_length=0)
        if mode == "perfect":
            ps._run_perfect(tel, seq, cfg, raw, max_workers=w)
        else:
            ps._run_imperfect(tel, seq, cfg, raw, mode, max_workers=w)
