#!/usr/bin/env python3
"""Unit timecourse figure + table from existing exact-search results (no new searches).

For each chosen unit: % of each sample's telomere ends with >= 1 exact copy (top panel) and
with >= 2 exact tandem copies (bottom panel). Reads <sample>_<unit>_perfect_search.txt (telomere
total) and <sample>_<unit>_perfect_search_matches.fasta (copies per read) from SEARCH_DIR.

Usage: plot_units_timecourse.py SEARCH_DIR SAMPLES_FILE STRAIN OUT_PREFIX UNIT=LABEL [UNIT=LABEL ...]
"""
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D, SF, STRAIN, OUT = sys.argv[1:5]
UNITS = [a.split("=", 1) for a in sys.argv[5:]]
SAMPLES = [l.strip() for l in open(SF) if l.strip()]
LAB = [s.replace(f"{STRAIN}_", "").replace("_with_selection", "").replace("survivor_", "")
       .replace("_repeat2", " rep2").replace("_repeat", " rep").replace("_TeloTag", " TeloTag") for s in SAMPLES]
COLS = ["#2a7f7f", "#c0502e", "#6a5acd", "#d4a017", "#3b7dd8", "#8c564b", "#999999"]

tab = open(f"{OUT}_counts.tsv", "w")
tab.write("sample\tunit\tlabel\ttelomere_ends\treads_ge1_exact\treads_ge2_tandem\tpct_ge1\tpct_ge2\tmax_tandem_copies\n")
fig, (a, b) = plt.subplots(2, 1, figsize=(max(13, 0.8 * len(SAMPLES)), 7.5))
w = 0.8 / len(UNITS)
for k, (u, label) in enumerate(UNITS):
    p1, p2, c1, c2 = [], [], [], []
    for s in SAMPLES:
        tel = int(re.search(r"Telomeric reads\s*:\s*(\d+)", open(f"{D}/{s}_{u}_perfect_search.txt").read()).group(1))
        c = [int(x) for x in re.findall(r"copies=(\d+)", open(f"{D}/{s}_{u}_perfect_search_matches.fasta").read())]
        n1, n2 = len(c), sum(x >= 2 for x in c)
        p1.append(100 * n1 / tel)
        p2.append(100 * n2 / tel)
        c1.append(n1)
        c2.append(n2)
        tab.write(f"{s}\t{u}\t{label}\t{tel}\t{n1}\t{n2}\t{100 * n1 / tel:.4f}\t{100 * n2 / tel:.4f}\t{max(c) if c else 0}\n")
    xs = [i + (k - (len(UNITS) - 1) / 2) * w for i in range(len(SAMPLES))]
    a.bar(xs, p1, w, color=COLS[k % len(COLS)], label=label)
    b.bar(xs, p2, w, color=COLS[k % len(COLS)], label=label)
    for ax, ps, cs in ((a, p1, c1), (b, p2, c2)):
        for x, p, n in zip(xs, ps, cs):
            if n:
                ax.text(x, p * 1.15, str(n), ha="center", va="bottom", fontsize=5.5, rotation=90)
tab.close()
for ax, t in ((a, "telomere ends with >= 1 exact copy"), (b, "telomere ends with >= 2 exact tandem copies (circle signature)")):
    ax.set_xticks(range(len(SAMPLES)))
    ax.set_xticklabels(LAB, fontsize=7, rotation=35, ha="right")
    ax.set_ylabel("% of all telomere ends (log)")
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 300)
    ax.set_title(f"{STRAIN}: {t} (no length gate)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    days = [re.search(r"day(\d+)", s).group(1) if re.search(r"day(\d+)", s) else "S" for s in SAMPLES]
    for i in range(1, len(days)):
        if days[i] != days[i - 1]:
            ax.axvline(i - 0.5, color="k", lw=0.5, ls=":")
a.legend(frameon=False, fontsize=8, ncol=len(UNITS))
fig.text(0.01, 0.005, "Exact matching only (nanopore errors undercount long units). Numbers above bars = reads. Log scale; missing bar = 0 reads. Dotted lines separate days.", fontsize=7)
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(f"{OUT}.png", dpi=200)
print(f"{OUT}.png\n{OUT}_counts.tsv")
