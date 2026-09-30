#!/usr/bin/env python3
"""Heatmap of every circle candidate population mode flagged (all families >= 90 bp) in one strain,
searched (exact matches, no length gate) in all of that strain's samples.
Usage: plot_all_candidates.py OUTDIR SAMPLES_FILE STRAIN

Rows = candidate units, grouped into "same circle" groups (>= 97% rotation-aware containment,
i.e. dimers / variants of one another), labelled by type. Columns = samples.
Top panel: telomere ends with >= 1 exact copy; bottom: >= 2 exact tandem copies.
Cell text = number of reads; colour = % of that sample's telomere ends (log-ish scale).
A black box marks the sample in which the program called that candidate."""
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

F, SF, STRAIN = sys.argv[1], sys.argv[2], sys.argv[3]
SAMPLES = [l.strip() for l in open(SF) if l.strip()]
SHORT = [re.sub(r"^\d+_", "", s).replace(f"{STRAIN}_", "").replace("_with_selection", "").replace("survivor_", "")
         .replace("_repeat", " rep") for s in SAMPLES]
units = [l.rstrip("\n").split("\t") for l in open(f"{F}/units.tsv")]
info = {r.split("\t")[0]: r.rstrip("\n").split("\t") for r in list(open(f"{F}/units_info.tsv"))[1:]}


def semi(a, B):
    prev = [0] * (len(B) + 1)
    for ca in a:
        cur = [prev[0] + 1] + [0] * len(B)
        for j, cb in enumerate(B, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
        prev = cur
    return min(prev)


def same(a, b):
    x, y = (a, b) if len(a) <= len(b) else (b, a)
    return 1 - semi(x, y + y) / len(x) >= 0.97


# single-linkage grouping of units into "same circle" groups
n = len(units)
parent = list(range(n))
def find(i):
    while parent[i] != i:
        parent[i] = parent[parent[i]]
        i = parent[i]
    return i
for i in range(n):
    for j in range(i + 1, n):
        if same(units[i][3], units[j][3]):
            parent[find(i)] = find(j)
groups = {}
for i in range(n):
    groups.setdefault(find(i), []).append(i)
# order: TG telomeric groups first (by first-called sample), then simple repeats
def gkey(g):
    simple = info[units[g[0]][0]][3].startswith("simple")
    return (simple, min(SAMPLES.index(units[i][1]) for i in g), min(len(units[i][3]) for i in g))
order = [i for g in sorted(groups.values(), key=gkey) for i in sorted(g, key=lambda i: (SAMPLES.index(units[i][1]), len(units[i][3])))]
gid = {}
for k, g in enumerate(sorted(groups.values(), key=gkey), 1):
    for i in g:
        gid[i] = k

# counts
tel = {}
one = [[0] * len(SAMPLES) for _ in range(n)]
two = [[0] * len(SAMPLES) for _ in range(n)]
for c, s in enumerate(SAMPLES):
    for r, u in enumerate(units):
        txt = open(f"{F}/{s}_{u[0]}_perfect_search.txt").read()
        tel[s] = int(re.search(r"Telomeric reads\s*:\s*(\d+)", txt).group(1))
        cp = [int(x) for x in re.findall(r"copies=(\d+)", open(f"{F}/{s}_{u[0]}_perfect_search_matches.fasta").read())]
        one[r][c] = len(cp)
        two[r][c] = sum(x >= 2 for x in cp)

with open(f"{F}/{STRAIN}_all{n}_counts.tsv", "w") as out:
    out.write("group\tunit\tlength\tkind\tcalled_in\t" + "\t".join(f"{x}_ge1\t{x}_ge2" for x in SHORT) + "\n")
    for i in order:
        u = units[i]
        out.write(f"{gid[i]}\t{u[0]}\t{len(u[3])}\t{info[u[0]][3]}\t{u[1]}\t"
                  + "\t".join(f"{one[i][c]}\t{two[i][c]}" for c in range(len(SAMPLES))) + "\n")

fig, axes = plt.subplots(2, 1, figsize=(15, 2 + 0.36 * n * 2))
for ax, M, title in ((axes[0], one, "telomere ends with >= 1 exact copy"),
                     (axes[1], two, "telomere ends with >= 2 exact tandem copies (circle signature)")):
    pct = [[100 * M[i][c] / tel[s] if M[i][c] else float("nan") for c, s in enumerate(SAMPLES)] for i in order]
    im = ax.imshow(pct, aspect="auto", cmap="viridis", norm=LogNorm(vmin=0.002, vmax=1))
    for r, i in enumerate(order):
        for c in range(len(SAMPLES)):
            v = M[i][c]
            if v:
                ax.text(c, r, str(v), ha="center", va="center", fontsize=6,
                        color="white" if pct[r][c] < 0.05 else "black")
        c0 = SAMPLES.index(units[i][1])
        ax.add_patch(Rectangle((c0 - 0.5, r - 0.5), 1, 1, fill=False, lw=1.2, ec="black"))
    labels = []
    for i in order:
        u = units[i]
        kind = info[u[0]][3]
        tag = "SIMPLE " + kind.split("(")[1].rstrip(")") if kind.startswith("simple") else "TG"
        labels.append(f"G{gid[i]:<2d} {u[0]:<18s} {tag}")
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=6.5, family="monospace")
    for t, i in zip(ax.get_yticklabels(), order):
        if info[units[i][0]][3].startswith("simple"):
            t.set_color("#a33")
    ax.set_xticks(range(len(SAMPLES)))
    ax.set_xticklabels([f"{x}\n({tel[s]:,})" for x, s in zip(SHORT, SAMPLES)], fontsize=7)
    surv = [k for k, s in enumerate(SAMPLES) if "survivor" in s]
    if surv:
        ax.axvline(surv[0] - 0.5, color="white", lw=2)
    prev = None
    for r, i in enumerate(order):
        if prev is not None and gid[i] != prev:
            ax.axhline(r - 0.5, color="white", lw=0.8)
        prev = gid[i]
    ax.set_title(f"{STRAIN}, all {n} circle candidates the program flagged: {title}. "
                 "Number = reads; colour = % of telomere ends; box = sample where the program called it", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.015, pad=0.01, label="% of telomere ends")
fig.text(0.01, 0.003, "Rows grouped (G#) when units match each other >= 97% (same circle / dimer / variant). "
         "Red labels = simple repeats (CTTT/AGGG/GTT...), not telomeric circles. Exact matching only: long units "
         "are undercounted by nanopore errors. Column header numbers = telomere ends searched (no length gate).",
         fontsize=6.5)
fig.tight_layout(rect=(0, 0.01, 1, 1))
fig.savefig(f"{F}/{STRAIN}_all{n}_circles_heatmap.png", dpi=200)
print(f"{n} units in {len(groups)} groups -> {F}/{STRAIN}_all{n}_circles_heatmap.png")
for k, g in enumerate(sorted(groups.values(), key=gkey), 1):
    print(f"G{k}: " + ", ".join(units[i][0] for i in g) + f"   [{info[units[g[0]][0]][3]}]")
