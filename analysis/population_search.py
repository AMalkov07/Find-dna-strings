"""
Population mode: search for a user-specified pattern in telomeric reads.

Three modes (selected via -psm / --population_search_mode):
  perfect            — exact tandem copy count per read; outputs summary,
                       copy-count histogram, and a FASTA of matching reads.
  alignment          — full alignment analysis on the top N reads by copy count.
  template_switching — full template-switching analysis on the top N reads.
"""

import copy
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from os.path import splitext
from typing import Dict, List, Optional, Tuple

from utils.data_structures import Config, TelomereSequence

# ── core helpers ──────────────────────────────────────────────────────────────

def _perfect_scan(sequence: str, pattern: str) -> Tuple[List[int], int]:
    """
    Find all non-overlapping exact occurrences of pattern in sequence.
    Returns (positions, max_consecutive_tandem_copies).
    """
    n = len(pattern)
    positions: List[int] = []
    i = 0
    while i <= len(sequence) - n:
        if sequence[i:i + n] == pattern:
            positions.append(i)
            i += n
        else:
            i += 1
    if not positions:
        return [], 0
    max_run = cur_run = 1
    for j in range(1, len(positions)):
        if positions[j] == positions[j - 1] + n:
            cur_run += 1
            if cur_run > max_run:
                max_run = cur_run
        else:
            cur_run = 1
    return positions, max_run


def _scan_one_read(args: Tuple) -> Optional[Tuple[str, List[int], int]]:
    """Module-level worker for ProcessPoolExecutor — scan one read for pattern."""
    read_id, sequence, pattern = args
    positions, max_tandem = _perfect_scan(sequence, pattern)
    if positions:
        return (read_id, positions, max_tandem)
    return None


def _rank_reads(sequences: Dict[str, str], pattern: str,
                max_workers: Optional[int] = None) -> List[Tuple[str, List[int], int]]:
    """
    Scan all reads for the pattern in parallel.
    Returns [(read_id, positions, max_tandem), ...] sorted by max_tandem descending.
    Only reads with at least one occurrence are included.
    """
    args = [(read_id, seq, pattern) for read_id, seq in sequences.items()]
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        raw = list(executor.map(_scan_one_read, args))
    results = [r for r in raw if r is not None]
    results.sort(key=lambda x: x[2], reverse=True)
    return results


def _write_histogram(f, copy_counts: List[int], max_bar: int = 40) -> None:
    dist = Counter(copy_counts)
    if not dist:
        return
    max_count = max(dist.values())
    f.write(f"\n  {'Copies':>7}  {'Reads':>6}  Bar\n")
    f.write("  " + "-" * 55 + "\n")
    for copies in sorted(dist.keys()):
        count = dist[copies]
        bar = '\u2588' * int(count / max_count * max_bar)
        f.write(f"  {copies:>7}  {count:>6}  {bar}\n")


# ── perfect match mode ────────────────────────────────────────────────────────

def _run_perfect(sequences: Dict[str, str], pattern: str,
                 config: Config, total_raw: int,
                 max_workers: Optional[int] = None) -> None:
    base, _ = splitext(config.output_file)
    out_path   = f"{base}_search.txt"
    fasta_path = f"{base}_search_matches.fasta"

    ranked = _rank_reads(sequences, pattern, max_workers=max_workers)
    n_tel  = len(sequences)
    copy_counts = [t for _, _, t in ranked]

    with open(out_path, 'w') as f:
        f.write("POPULATION PATTERN SEARCH — PERFECT MATCH\n")
        f.write("=" * 70 + "\n")
        f.write(f"  Pattern         : {pattern}\n")
        f.write(f"  Pattern length  : {len(pattern)} bp\n")
        f.write(f"  Total raw reads : {total_raw}\n")
        f.write(f"  Telomeric reads : {n_tel}\n")
        f.write("=" * 70 + "\n\n")

        f.write("MATCH SUMMARY\n")
        f.write("-" * 70 + "\n")
        for threshold in [1, 2, 3, 4, 5]:
            n = sum(1 for c in copy_counts if c >= threshold)
            label = "copy " if threshold == 1 else "copies"
            f.write(f"  Reads with \u2265{threshold} {label}: {n:>6} / {n_tel}  ({n/n_tel:.1%})\n")
        if copy_counts:
            f.write(f"\n  Max tandem copies  : {max(copy_counts)}\n")
            f.write(f"  Mean tandem copies : {sum(copy_counts)/len(copy_counts):.2f}"
                    f"  (positive reads only)\n")

        f.write("\nTANDEM COPY COUNT DISTRIBUTION\n")
        f.write("-" * 70 + "\n")
        _write_histogram(f, copy_counts)
        f.write(f"\nMatching reads written to: {fasta_path}\n")

    with open(fasta_path, 'w') as f:
        for read_id, positions, max_tandem in ranked:
            f.write(f">{read_id}  copies={max_tandem}  occurrences={len(positions)}\n"
                    f"{sequences[read_id]}\n")

    print(f"  Perfect match : {len(ranked)}/{n_tel} reads contain the pattern")
    print(f"  Results       : {out_path}")
    print(f"  Matches FASTA : {fasta_path}")


# ── alignment / template-switching mode ──────────────────────────────────────

def _run_imperfect(sequences: Dict[str, str], pattern: str,
                   config: Config, total_raw: int,
                   mode: str, top_n: int = 20,
                   max_workers: Optional[int] = None) -> None:
    from analysis.alignment_strategy import AlignmentStrategy
    from analysis.template_switching_strategy_2 import TemplateSwitchingStrategy
    from data_io.template_switching_exporters import TemplateSwitchingPrint
    from data_io.alignment_exporters import AlignmentPrint

    base, _ = splitext(config.output_file)
    out_path = f"{base}_search.txt"

    ranked = _rank_reads(sequences, pattern, max_workers=max_workers)
    n_tel  = len(sequences)
    top    = ranked[:top_n]

    telomere_objects = [
        TelomereSequence(
            survivor_name=read_id,
            sequence=sequences[read_id],
            chromosome_end_id=read_id,
            analysis=None,
        )
        for read_id, _, _ in top
    ]

    search_config = copy.copy(config)
    search_config.output_file = out_path
    search_config.graph_output = None   # no graphs for population search

    # Write header (file opened with 'w' so it's fresh)
    with open(out_path, 'w') as f:
        label = mode.upper().replace('_', ' ')
        f.write(f"POPULATION PATTERN SEARCH — {label}\n")
        f.write("=" * 70 + "\n")
        f.write(f"  Pattern         : {pattern}\n")
        f.write(f"  Pattern length  : {len(pattern)} bp\n")
        f.write(f"  Total raw reads : {total_raw}\n")
        f.write(f"  Telomeric reads : {n_tel}\n")
        f.write(f"  Reads with \u22651 perfect copy : {len(ranked)} ({len(ranked)/n_tel:.1%})\n")
        f.write(f"  Showing top {len(top)} reads by tandem copy count\n")
        f.write("=" * 70 + "\n\n")

    # Run strategy then append detailed per-read output
    if mode == 'template_switching':
        analyzer = TemplateSwitchingStrategy(telomere_objects, pattern, search_config)
        analyzer.execute()
        TemplateSwitchingPrint(telomere_objects, search_config, pattern).print_analysis(file_mode='a')
    else:
        analyzer = AlignmentStrategy(telomere_objects, pattern, search_config)
        analyzer.execute()
        # Build per-read copy counts (perfect + imperfect alignments) for histogram
        copy_counts = []
        for t in telomere_objects:
            if t.analysis is not None:
                n = len(t.analysis.perfect_alignments) + len(t.analysis.imperfect_alignments)
                if n > 0:
                    copy_counts.append(n)
        with open(out_path, 'a') as f:
            f.write("TANDEM COPY COUNT DISTRIBUTION (perfect + imperfect alignments)\n")
            f.write("-" * 70 + "\n")
            _write_histogram(f, copy_counts)
            f.write("\n")
        AlignmentPrint(telomere_objects, search_config, pattern).print_analysis(file_mode='a')

    print(f"  {mode} search : top {len(top)} of {len(ranked)} positive reads analyzed")
    print(f"  Results       : {out_path}")


# ── entry point ───────────────────────────────────────────────────────────────

def population_search_execute(config: Config) -> None:
    """Search for a specific pattern in population telomeric reads."""
    from data_io.fasta_reader import FastaReader
    from analysis.pattern_finder_lcp import _resolve_workers

    if not config.pattern:
        raise ValueError("population search requires a pattern via -p / --pattern")

    max_workers = _resolve_workers(config.workers)
    reader = FastaReader(config.fasta_file_path, max_ends=0)
    sequences, total_raw = reader.parse_fasta_population(
        threshold=config.telomere_threshold, max_workers=max_workers
    )

    if not sequences:
        raise ValueError("No telomeric sequences found in FASTA file")

    print(f"  Telomere extraction: {len(sequences)} telomeric reads from {total_raw} total reads")

    mode = config.population_search_mode
    if mode == 'perfect':
        _run_perfect(sequences, config.pattern, config, total_raw, max_workers=max_workers)
    else:
        _run_imperfect(sequences, config.pattern, config, total_raw, mode, max_workers=max_workers)
