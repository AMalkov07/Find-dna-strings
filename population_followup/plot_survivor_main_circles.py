#!/usr/bin/env python3
"""One bar chart, all survivors: each survivor's own main circle.

For each survivor, the main circle = the survivor's own TG family (>= 100 bp) with the most tandem ("circularized") reads.
Bar = % of that survivor's analysed telomeres (>= 200 bp) carrying the circle at >= 98% identity (light) and
with >= 2 tandem copies (dark), from classify_family_reads.py (exact matching undercounts some units);
label = circle length (bp) and read counts. Colour = circle group (units matching >= 97%, from the
heatmap run), so survivors that share a circle share a colour.

Reads OUTDIR/<strain>_all<N>_counts.tsv and OUTDIR/units_info.tsv (from plot_all_candidates.py /
build_units.py), plus each sample's population-mode outputs one level above OUTDIR:
<parent>/<sample>/<sample>_output_pattern.txt and _output_family_reads.txt.
Writes OUTDIR/<strain>_main_circle_per_survivor.{png,tsv}.

Usage: plot_survivor_main_circles.py OUTDIR SAMPLES_FILE STRAIN
"""
import csv
import glob
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

F, SF, STRAIN = sys.argv[1:4]
SAMPLES = [l.strip() for l in open(SF) if l.strip()]
short = lambda s: re.sub(r"^\d+_", "", s).replace("survivor_", "")
counts = list(csv.DictReader(open(glob.glob(f"{F}/{STRAIN}_all*_counts.tsv")[0]), delimiter="\t"))
info = {r["unit"]: r for r in csv.DictReader(open(f"{F}/units_info.tsv"), delimiter="\t")}
col_of = lambda s, k: next(c for c in counts[0] if c.endswith(f"{short(s)}_{k}"))

W = F.rsplit("/", 1)[0]
rows = []
for s in SAMPLES:
    pat = open(f"{W}/{s}/{s}_output_pattern.txt").read()
    tel = int(re.search(r"Telomeric reads\s*:\s*(\d+)", pat).group(1))          # telomeres >= 200 bp analysed
    fr = f"{W}/{s}/{s}_output_family_reads.txt"
    fams = re.findall(r"FAMILY (\d+)\s+\|\s+representative (\d+) bp\s+\|\s+circularized: (\d+)\s+\|\s+multi-copy: (\d+)"
                      r"\s+\|\s+single: (\d+)", open(fr).read()) if __import__("os").path.exists(fr) else []
    unit_of = {(r["sample"], r["family"]): u for u, r in info.items()}
    tg = [(int(f), int(L), int(c), int(m) + int(sg)) for f, L, c, m, sg in fams
          if (s, f"family{f}") in unit_of and not info[unit_of[(s, f"family{f}")]]["kind"].startswith("simple")]
    if tg:
        f, L, circ, other = max(tg, key=lambda x: (x[2], x[3]))
        u = unit_of[(s, f"family{f}")]
        g = next(r["group"] for r in counts if r["unit"] == u)
        rows.append((s, L, g, circ + other, circ, tel, u))
    else:
        called = [u for u, r in info.items() if r["sample"] == s]
        rows.append((s, None, "simple repeat only" if called else "no candidate", 0, 0, tel, ""))

groups = sorted({r[2] for r in rows if r[1]}, key=lambda g: int(g))
cmap = plt.get_cmap("tab20")
gcol = {g: cmap(i % 20) for i, g in enumerate(groups)}
shared = {g for g in groups if sum(1 for r in rows if r[2] == g) > 1}

with open(f"{F}/{STRAIN}_main_circle_per_survivor.tsv", "w") as out:
    out.write("sample\tmain_circle_unit\tlength_bp\tgroup\ttelomere_ends\treads_with_circle_98pct\treads_tandem_98pct\tpct_with_circle\tpct_tandem\tshared_with\n")
    for s, L, g, n1, n2, tel, u in rows:
        sh = ", ".join(short(r[0]) for r in rows if r[2] == g and r[0] != s) if L else ""
        out.write(f"{s}\t{u}\t{L or ''}\t{g}\t{tel or ''}\t{n1}\t{n2}\t"
                  f"{100 * n1 / tel if tel else 0:.4f}\t{100 * n2 / tel if tel else 0:.4f}\t{sh}\n")

fig, ax = plt.subplots(figsize=(max(14, 0.62 * len(rows)), 6.2))
for i, (s, L, g, n1, n2, tel, u) in enumerate(rows):
    if not L:
        ax.text(i, 1, f"{g}  ({tel:,} telomeres)", rotation=90, ha="center", va="bottom", fontsize=7, color="#888")
        continue
    p1, p2 = 100 * n1 / tel, 100 * n2 / tel
    c = gcol[g]
    if p1:
        ax.bar(i, p1, 0.8, color=c, alpha=0.45)
    if p2:
        ax.bar(i, p2, 0.8, color=c)
    top = max(p1, p2)
    ax.text(i, top + 1.5, f"{L} bp  ({n1:,} / {n2:,})", rotation=90, ha="center", va="bottom", fontsize=7,
            fontweight="bold" if g in shared else "normal")
ax.set_ylim(0, 145)
ax.set_yticks(range(0, 101, 20))
ax.set_xticks(range(len(rows)))
ax.set_xticklabels([short(r[0]) for r in rows], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("% of the survivor's telomeres >= 200 bp")
ax.set_title(f"{STRAIN}: each survivor's main telomeric circle. Light = reads carrying it, dark = reads with >= 2 tandem "
             "copies (>= 98% identity). Label = circle length (reads carrying / tandem). Same colour = same circle (>= 97% identity)",
             fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(handles=[Patch(color=gcol[g], label=f"G{g}: " + ", ".join(short(r[0]) for r in rows if r[2] == g))
                   for g in groups if g in shared],
          title="circles shared by several survivors", fontsize=7, title_fontsize=8, frameon=False,
          loc="upper left", ncol=2)
fig.text(0.01, 0.005, "Counts from classify_family_reads.py (>= 98% identity, full-length copies) over each survivor's telomeres >= 200 bp. "
         "Main circle = the survivor's own TG family (>= 100 bp) with the most tandem reads. Bold labels = circle shared with another survivor.", fontsize=7)
fig.tight_layout(rect=(0, 0.02, 1, 1))
fig.savefig(f"{F}/{STRAIN}_main_circle_per_survivor.png", dpi=200)
print(f"{F}/{STRAIN}_main_circle_per_survivor.png")
