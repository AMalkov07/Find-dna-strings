#!/usr/bin/env python3
"""
Post-process a population-mode *_output_pattern.txt file and mark patterns that
are variations of a higher-ranked pattern.

A pattern P is a "variation" of a higher-ranked pattern Q when, treating both as
circular sequences, the shorter of the two fits inside the other (rotation-aware,
allowing sequencing-error mismatches/indels) at >= IDENTITY_THRESHOLD identity.
Each pattern is attributed to the HIGHEST-RANKED pattern it is a direct
variation of. This covers both cases the user asked for:
  - P is a shorter variation that fits inside a higher-ranked pattern, and
  - P is a longer variation that a higher-ranked pattern fits inside of.

Pure standard library (no edlib/parasail/biopython needed). Patterns are few and
short, so the O(n*m) fuzzy-substring DP is plenty fast.

Usage:
    python annotate_pattern_variations.py FILE [FILE ...]
    python annotate_pattern_variations.py *_output_pattern.txt
    python annotate_pattern_variations.py --threshold 0.92 FILE
    python annotate_pattern_variations.py --in-place FILE        # overwrite

By default writes "<name>_annotated.txt" next to each input file.
"""

import argparse
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

# Shared containment logic (single source of truth, also used by the finder).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysis.pattern_variations import (  # noqa: E402
    DEFAULT_THRESHOLD, variant_identity,
)


# ── parsing ───────────────────────────────────────────────────────────────────

# Matches the per-rank blocks in the FULL PATTERN SEQUENCES section:
#   "  Rank 3  (278 bp  reads=20/34461  ...)\n    GGGTGT...\n"
_FULL_BLOCK_RE = re.compile(
    r"^\s*Rank\s+(\d+)\s+\(.*?\n\s+([ACGTNacgtn]+)\s*$",
    re.MULTILINE,
)


def parse_patterns(text: str) -> List[Tuple[int, str]]:
    """Return [(rank, sequence), ...] from the FULL PATTERN SEQUENCES section."""
    if "FULL PATTERN SEQUENCES" in text:
        section = text.split("FULL PATTERN SEQUENCES", 1)[1]
    else:
        section = text
    pats = [(int(r), seq.upper()) for r, seq in _FULL_BLOCK_RE.findall(section)]
    pats.sort(key=lambda x: x[0])
    return pats


# ── variation assignment ──────────────────────────────────────────────────────

def assign_variations(pats: List[Tuple[int, str]],
                      threshold: float) -> Dict[int, Optional[Tuple[int, float]]]:
    """Map each rank -> (higher_rank, identity) it is a variation of, or None.

    Patterns are assumed already ordered best-first. Each pattern is matched to
    the highest-ranked earlier pattern that meets the identity threshold.
    """
    result: Dict[int, Optional[Tuple[int, float]]] = {}
    for idx, (rank, seq) in enumerate(pats):
        match: Optional[Tuple[int, float]] = None
        for jdx in range(idx):  # earlier == higher-ranked
            hrank, hseq = pats[jdx]
            ident = variant_identity(seq, hseq)
            if ident >= threshold:
                match = (hrank, ident)
                break  # first hit is the highest-ranked one
        result[rank] = match
    return result


def _note(variation: Optional[Tuple[int, float]]) -> str:
    if variation is None:
        return ""
    hrank, ident = variation
    return f"   <- variation of Rank {hrank} ({ident:.1%} id)"


def build_summary(pats: List[Tuple[int, str]],
                  variations: Dict[int, Optional[Tuple[int, float]]],
                  threshold: float) -> str:
    """A human-readable summary block grouping patterns into families."""
    lines = []
    lines.append("=" * 90)
    lines.append("PATTERN VARIATION ANALYSIS")
    lines.append("=" * 90)
    lines.append(f"  Rotation-aware containment, >= {threshold:.0%} identity "
                 f"(allows mismatches/indels).")
    lines.append("  A pattern is a 'variation' of the highest-ranked pattern it "
                 "fits inside of (or that fits inside it).")
    lines.append("")
    lines.append(f"  {'Rank':>4}  {'Length':>6}  Variation")
    lines.append("  " + "-" * 60)
    length_by_rank = {rank: len(seq) for rank, seq in pats}
    for rank, seq in pats:
        v = variations[rank]
        if v is None:
            desc = "primary"
        else:
            hrank, ident = v
            desc = f"variation of Rank {hrank} ({ident:.1%} id)"
        lines.append(f"  {rank:>4}  {len(seq):>6}  {desc}")

    # Families: primary -> [variant ranks]
    families: Dict[int, List[int]] = {}
    for rank, seq in pats:
        v = variations[rank]
        if v is None:
            families.setdefault(rank, families.get(rank, []))
        else:
            families.setdefault(v[0], []).append(rank)
    lines.append("")
    lines.append("  Families (primary pattern -> its variations):")
    primaries = [rank for rank, _ in pats if variations[rank] is None]
    for prank in primaries:
        variants = sorted(families.get(prank, []))
        if variants:
            vtxt = ", ".join(str(r) for r in variants)
            lines.append(f"    Rank {prank} ({length_by_rank[prank]} bp): variations -> {vtxt}")
        else:
            lines.append(f"    Rank {prank} ({length_by_rank[prank]} bp): no variations")
    lines.append("=" * 90)
    lines.append("")
    return "\n".join(lines)


# ── annotation of the original text ───────────────────────────────────────────

def annotate_text(text: str, threshold: float) -> Tuple[str, int]:
    """Return (annotated_text, num_variations_found)."""
    pats = parse_patterns(text)
    if not pats:
        raise ValueError("no 'FULL PATTERN SEQUENCES' patterns found in file")
    variations = assign_variations(pats, threshold)
    n_var = sum(1 for v in variations.values() if v is not None)
    summary = build_summary(pats, variations, threshold)

    lines = text.splitlines(keepends=True)
    out: List[str] = []
    section = None  # "TOP" | "FULL" | None
    top_header_seen = False

    # rank -> annotation note
    notes = {rank: _note(variations[rank]) for rank, _ in pats}

    for line in lines:
        stripped = line.rstrip("\n")

        # Track sections.
        if "TOP CANDIDATE PATTERNS" in line:
            section = "TOP"
            top_header_seen = False
            out.append(line)
            continue
        if "FULL PATTERN SEQUENCES" in line:
            # Insert the summary just before the full-sequence section.
            out.append(summary)
            section = "FULL"
            out.append(line)
            continue
        if "CIRCLE DETECTION RESULT" in line:
            section = None
            out.append(line)
            continue

        if section == "TOP":
            # Data rows start with an integer rank; skip the column header row.
            m = re.match(r"^\s*(\d+)\s+\d+\s", line)
            if m:
                rank = int(m.group(1))
                if rank in notes and notes[rank]:
                    out.append(stripped + notes[rank] + "\n")
                    continue
            out.append(line)
            continue

        if section == "FULL":
            m = re.match(r"^\s*Rank\s+(\d+)\b", line)
            if m:
                rank = int(m.group(1))
                if rank in notes and notes[rank]:
                    out.append(stripped + notes[rank] + "\n")
                    continue
            out.append(line)
            continue

        out.append(line)

    return "".join(out), n_var


# ── cli ───────────────────────────────────────────────────────────────────────

def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="*_output_pattern.txt file(s) to annotate")
    ap.add_argument("-t", "--threshold", type=float, default=DEFAULT_THRESHOLD,
                    help=f"identity threshold for a variation match (default {DEFAULT_THRESHOLD})")
    ap.add_argument("--in-place", action="store_true",
                    help="overwrite the input file instead of writing *_annotated.txt")
    ap.add_argument("--suffix", default="_annotated",
                    help="suffix for output files when not --in-place (default: _annotated)")
    args = ap.parse_args(argv)

    rc = 0
    for path in args.files:
        try:
            with open(path) as f:
                text = f.read()
            annotated, n_var = annotate_text(text, args.threshold)
        except (OSError, ValueError) as e:
            print(f"  SKIP {path}: {e}", file=sys.stderr)
            rc = 1
            continue

        if args.in_place:
            out_path = path
        elif path.endswith(".txt"):
            out_path = path[:-4] + args.suffix + ".txt"
        else:
            out_path = path + args.suffix

        with open(out_path, "w") as f:
            f.write(annotated)
        print(f"  {path}  ->  {out_path}   ({n_var} variation(s) marked)")
    return rc


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
