import math
import os
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor
from typing import List, Tuple, Dict, Set, Optional
from utils.data_structures import TelomereSequence, Config

from os.path import splitext
from _io import TextIOWrapper


def _resolve_workers(workers: Optional[int]) -> Optional[int]:
    """Return explicit workers count, or read NSLOTS (SGE), or None (let ProcessPoolExecutor decide)."""
    if workers is not None:
        return workers
    nslots = os.environ.get('NSLOTS')
    return int(nslots) if nslots else None


def _analyze_single_sequence(args: Tuple) -> Tuple[str, Dict]:
    """
    Module-level worker for ProcessPoolExecutor — analyze one sequence for repeated patterns.
    Must be module-level (not a method) so multiprocessing can pickle it.
    """
    seq_id, sequence, min_length, max_length = args
    n = len(sequence)

    # ── Suffix array ──────────────────────────────────────────────────────────
    suffix_array = sorted(range(n), key=lambda i: sequence[i:])

    # ── LCP array ─────────────────────────────────────────────────────────────
    lcp_array = []
    for i in range(n - 1):
        s1, s2 = suffix_array[i], suffix_array[i + 1]
        common = 0
        limit = min(n - s1, n - s2)
        for j in range(limit):
            if sequence[s1 + j] == sequence[s2 + j]:
                common += 1
            else:
                break
        lcp_array.append(common)

    # ── Raw repeated patterns via LCP ─────────────────────────────────────────
    raw_patterns: Dict[str, set] = defaultdict(set)
    for i, lcp_len in enumerate(lcp_array):
        if lcp_len >= min_length:
            pos1 = suffix_array[i]
            pos2 = suffix_array[i + 1]
            max_extract = min(lcp_len, max_length)
            for length in range(min_length, max_extract + 1):
                if pos1 + length <= n:
                    pattern = sequence[pos1:pos1 + length]
                    raw_patterns[pattern].add(pos1)
                    raw_patterns[pattern].add(pos2)

    # ── Consolidate rotational patterns ───────────────────────────────────────
    canonical_groups: Dict[str, list] = defaultdict(list)
    for pattern, position_set in raw_patterns.items():
        p_len = len(pattern)
        doubled = pattern + pattern
        canonical = min(doubled[i:i + p_len] for i in range(p_len))
        canonical_groups[canonical].append(position_set)

    repeated_patterns: Dict[str, List[int]] = {}
    for canonical, pos_sets in canonical_groups.items():
        all_pos: set = set()
        for ps in pos_sets:
            all_pos.update(ps)
        repeated_patterns[canonical] = sorted(all_pos)

    # ── Build per-sequence analysis ───────────────────────────────────────────
    sequence_analysis: Dict = {}
    for canonical_pattern, all_positions in repeated_patterns.items():
        if len(all_positions) < 2:
            continue
        pattern_len = len(canonical_pattern)

        # Non-overlapping positions (greedy)
        non_overlapping: List[int] = []
        last_end = -1
        for pos in all_positions:  # already sorted
            if pos >= last_end:
                non_overlapping.append(pos)
                last_end = pos + pattern_len

        # Tandem arrays
        tandem_arrays: List[Tuple[int, int, int]] = []
        sp = non_overlapping
        i = 0
        while i < len(sp):
            cur = sp[i]
            count = 1
            for j in range(i + 1, len(sp)):
                if sp[j] == cur + pattern_len:
                    count += 1
                    cur = sp[j]
                else:
                    break
            if count >= 2:
                tandem_arrays.append((sp[i], cur + pattern_len, count))
                i += count
            else:
                i += 1

        sequence_analysis[canonical_pattern] = {
            'positions': all_positions,
            'non_overlapping_positions': non_overlapping,
            'total_occurrences': len(all_positions),
            'non_overlapping_occurrences': len(non_overlapping),
            'tandem_arrays': tandem_arrays,
            'tandem_array_count': len(tandem_arrays),
            'total_tandem_copies': sum(c for _, _, c in tandem_arrays),
            'max_tandem_copies': max((c for _, _, c in tandem_arrays), default=0),
            'pattern_length': pattern_len,
        }

    return seq_id, sequence_analysis

class ChromosomeEndRepeatFinder:
    def __init__(self, sequences: Dict[str, str], pattern_file: TextIOWrapper,
                 max_workers: Optional[int] = None):
        """
        Initialize with multiple DNA sequences (e.g., chromosome ends)
        sequences: Dict mapping sequence_id -> sequence_string
        """
        self.sequences = {seq_id: seq.upper().replace('N', '') for seq_id, seq in sequences.items()}
        self.sequence_ids = list(self.sequences.keys())
        self.max_workers = max_workers
        print(f"Loaded {len(self.sequences)} sequences")
        self.pattern_file = pattern_file
        if len(self.sequences) <= 50:
            for seq_id, seq in self.sequences.items():
                pattern_file.write(f"  {seq_id:<8}  {len(seq):>10,} bp\n")
        else:
            lengths = [len(seq) for seq in self.sequences.values()]
            pattern_file.write(f"  (listing omitted — {len(self.sequences)} reads)\n")
            pattern_file.write(f"  Length range : {min(lengths):,} – {max(lengths):,} bp\n")
            pattern_file.write(f"  Mean length  : {sum(lengths) // len(lengths):,} bp\n")
    
    def _build_suffix_array(self, sequence: str) -> List[int]:
        """Build suffix array for a sequence"""
        n = len(sequence)
        suffixes = [(sequence[i:], i) for i in range(n)]
        suffixes.sort(key=lambda x: x[0])
        return [suffix[1] for suffix in suffixes]
    
    def _build_lcp_array(self, sequence: str, suffix_array: List[int]) -> List[int]:
        """Build LCP array from suffix array"""
        n = len(sequence)
        if n <= 1:
            return []
        
        lcp = []
        for i in range(n - 1):
            suffix1_start = suffix_array[i]
            suffix2_start = suffix_array[i + 1]
            
            common_len = 0
            max_compare = min(n - suffix1_start, n - suffix2_start)
            
            for j in range(max_compare):
                if sequence[suffix1_start + j] == sequence[suffix2_start + j]:
                    common_len += 1
                else:
                    break
            lcp.append(common_len)
        
        return lcp
    
    def get_canonical_rotation(self, pattern: str) -> str:
        """
        Get the lexicographically smallest rotation of a pattern
        This ensures all circular permutations are treated as the same pattern
        """
        if not pattern:
            return pattern
        n = len(pattern)
        doubled = pattern + pattern
        return min(doubled[i:i+n] for i in range(n))
    
    def is_rotation(self, str1: str, str2: str) -> bool:
        """
        Fast check if str2 is a rotation of str1 using the doubling trick
        """
        if len(str1) != len(str2):
            return False
        if str1 == str2:
            return True
        return str2 in (str1 + str1)
    
    def find_all_pattern_occurrences_fast(self, sequence: str, pattern: str) -> List[int]:
        """
        Fast pattern matching using doubling trick to find all occurrences (including rotations)
        """
        pattern_len = len(pattern)
        doubled_pattern = pattern + pattern
        all_positions = []
        
        # Find all occurrences in sequence
        start = 0
        while start <= len(sequence) - pattern_len:
            pos = sequence.find(pattern, start)
            if pos == -1:
                break
            all_positions.append(pos)
            start = pos + 1
        
        # Also find rotations by checking if substrings match any rotation
        # This is still expensive, so we'll do it differently...
        return all_positions
    
    def consolidate_rotational_patterns(self, raw_patterns: Dict[str, Set[int]], 
                                      sequence: str) -> Dict[str, List[int]]:
        """
        Consolidate rotational patterns early - keep only the best representative
        and find ALL positions for that representative (including rotations)
        """
        # Group patterns by their canonical rotation
        canonical_groups = defaultdict(list)  # canonical -> [(pattern, positions), ...]

        for pattern, position_set in raw_patterns.items():
            canonical = self.get_canonical_rotation(pattern)
            canonical_groups[canonical].append((pattern, position_set))
        
        # For each canonical group, union all positions already found for every rotation.
        # The LCP analysis records every position where each rotation appears (≥2 times),
        # so the union is complete — no need to rescan the full sequence again.
        final_patterns = {}

        for canonical, pattern_list in canonical_groups.items():
            all_positions: Set[int] = set()
            for _, positions in pattern_list:
                all_positions.update(positions)
            final_patterns[canonical] = sorted(all_positions)

        return final_patterns
    
    def find_all_rotations_efficiently(self, sequence: str, pattern: str) -> Set[int]:
        """
        Efficiently find all positions where pattern or its rotations occur
        Uses the doubling trick concept but optimized for sequence scanning
        """
        pattern_len = len(pattern)
        all_positions = set()
        
        # Create doubled pattern for rotation checking
        doubled_pattern = pattern + pattern
        
        # Scan sequence with sliding window
        for i in range(len(sequence) - pattern_len + 1):
            candidate = sequence[i:i + pattern_len]
            
            # Check if candidate is a rotation of pattern using doubling trick
            if candidate in doubled_pattern:
                all_positions.add(i)
        
        return all_positions
    
    def find_repeated_patterns_in_sequence(self, sequence: str, min_length: int, 
                                         max_length: int) -> Dict[str, List[int]]:
        """
        Use suffix array + LCP to find all repeated patterns, with early rotation consolidation
        Returns: Dict[canonical_pattern -> List[all_positions_including_rotations]]
        """
        suffix_array = self._build_suffix_array(sequence)
        lcp_array = self._build_lcp_array(sequence, suffix_array)
        raw_patterns = defaultdict(set)  # Use set to avoid duplicate positions
        
        # Process LCP array to find repeated substrings
        for i, lcp_len in enumerate(lcp_array):
            if lcp_len >= min_length:
                pos1 = suffix_array[i]
                pos2 = suffix_array[i + 1]
                
                # Extract patterns of different lengths from this common prefix
                max_extract = min(lcp_len, max_length)
                
                for length in range(min_length, max_extract + 1, 1):  # Larger step to reduce work
                    if pos1 + length <= len(sequence):
                        pattern = sequence[pos1:pos1 + length]
                        # Don't canonicalize yet - we'll do it in consolidation step
                        raw_patterns[pattern].add(pos1)
                        raw_patterns[pattern].add(pos2)
        
        # Consolidate rotational patterns and find all their positions
        return self.consolidate_rotational_patterns(raw_patterns, sequence)
    
    def find_tandem_arrays_optimized(self, known_positions: List[int], pattern_len: int) -> List[Tuple[int, int, int]]:
        """
        Find tandem arrays from known positions (much faster!)
        Args:
            known_positions: List of positions where pattern occurs (already found)
            pattern_len: Length of the pattern
        Returns: List[(start_pos, end_pos, copy_count)]
        """
        if len(known_positions) < 2:
            return []
        
        tandem_arrays = []
        sorted_positions = sorted(known_positions)
        
        i = 0
        while i < len(sorted_positions):
            tandem_start = sorted_positions[i]
            copy_count = 1
            current_pos = tandem_start
            
            # Count consecutive copies
            for j in range(i + 1, len(sorted_positions)):
                expected_pos = current_pos + pattern_len
                if sorted_positions[j] == expected_pos:
                    copy_count += 1
                    current_pos = expected_pos
                else:
                    break
            
            # If we found multiple consecutive copies, record this tandem array
            if copy_count >= 2:
                tandem_end = current_pos + pattern_len
                tandem_arrays.append((tandem_start, tandem_end, copy_count))
                i += copy_count  # Skip the positions we just processed
            else:
                i += 1
        
        return tandem_arrays
    
    def find_tandem_arrays(self, sequence: str, pattern: str) -> List[Tuple[int, int, int]]:
        """
        Find tandem arrays (consecutive repeats) of a large pattern
        Returns: List[(start_pos, end_pos, copy_count)]
        
        NOTE: This is the old slow method, kept for compatibility
        """
        pattern_len = len(pattern)
        tandem_arrays = []
        
        # Find all occurrences of the pattern
        all_positions = []
        for start in range(len(sequence) - pattern_len + 1):
            if sequence[start:start + pattern_len] == pattern:
                all_positions.append(start)
        
        return self.find_tandem_arrays_optimized(all_positions, pattern_len)
    
    def get_non_overlapping_positions(self, positions: List[int], pattern_len: int) -> List[int]:
        """
        Filter positions to get only non-overlapping occurrences
        Uses greedy algorithm: take first position, skip all that overlap, take next, etc.
        """
        if not positions:
            return []
        
        sorted_positions = sorted(positions)
        non_overlapping = []
        last_end = -1
        
        for pos in sorted_positions:
            if pos >= last_end:  # No overlap with previous
                non_overlapping.append(pos)
                last_end = pos + pattern_len
        
        return non_overlapping
    
    def analyze_all_sequences(self, min_length: int, max_length: int) -> Dict[str, Dict[str, any]]:
        """
        Analyze all sequences to find repeated patterns and their tandem arrays.
        Runs in parallel using ProcessPoolExecutor — one worker per sequence.
        Returns: Dict[seq_id -> Dict[pattern -> analysis_data]]
        """
        args = [(sid, seq, min_length, max_length) for sid, seq in self.sequences.items()]

        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            results = list(executor.map(_analyze_single_sequence, args))

        all_sequence_data = {}
        for seq_id, sequence_analysis in results:
            seq_len = len(self.sequences[seq_id])
            self.pattern_file.write(
                f"  Analyzed {seq_id} ({seq_len:,} bp): {len(sequence_analysis)} patterns\n"
            )
            all_sequence_data[seq_id] = sequence_analysis

        return all_sequence_data
    
    def count_cross_sequence_pattern_occurrences(self, target_pattern: str, sequence: str) -> int:
        """
        Fast count of pattern occurrences (including rotations) in a sequence
        Uses the doubling trick for efficient rotation matching
        """
        pattern_len = len(target_pattern)
        doubled_target = target_pattern + target_pattern
        count = 0
        
        # Scan sequence looking for the target pattern or any of its rotations
        for i in range(len(sequence) - pattern_len + 1):
            candidate = sequence[i:i + pattern_len]
            if candidate in doubled_target:
                count += 1
        
        return count
    
    def get_all_rotations(self, pattern: str) -> List[str]:
        """Get all circular rotations of a pattern"""
        rotations = []
        for i in range(len(pattern)):
            rotation = pattern[i:] + pattern[:i]
            rotations.append(rotation)
        return rotations

    def get_minimal_period(self, pattern: str) -> int:
        """
        Find the minimal period of a string using the KMP failure function.
        Returns the length of the shortest substring Q such that pattern = Q^k for some k >= 1.
        If the pattern is primitive (no shorter period divides evenly), returns len(pattern).
        """
        n = len(pattern)
        if n == 0:
            return 0
        fail = [0] * n
        j = 0
        for i in range(1, n):
            while j > 0 and pattern[i] != pattern[j]:
                j = fail[j - 1]
            if pattern[i] == pattern[j]:
                j += 1
            fail[i] = j
        period = n - fail[n - 1]
        return period if n % period == 0 else n

    def _filter_tandem_duplicates(self, scored_patterns: List, qualified_patterns: Set[str],
                                  min_length: int) -> List:
        """
        Remove patterns that are exact tandem repetitions of a shorter primitive pattern
        already present in the qualified candidate set.

        Example: a 130 bp pattern that equals [65 bp][65 bp] is dropped when the 65 bp
        canonical pattern is already a qualified candidate.
        """
        filtered = []
        for pattern, stats, score in scored_patterns:
            period = self.get_minimal_period(pattern)
            if period < len(pattern) and period >= min_length:
                # Pattern is a tandem repeat of its primitive root
                primitive_root = pattern[:period]
                canonical_root = self.get_canonical_rotation(primitive_root)
                if canonical_root in qualified_patterns:
                    repeats = len(pattern) // period
                    self.pattern_file.write(
                        f"    [filter] Dropping {len(pattern)} bp pattern "
                        f"(= {period} bp primitive × {repeats})\n"
                    )
                    continue
            filtered.append((pattern, stats, score))
        return filtered
    
    def merge_overlapping_arrays(self, arrays: List[Tuple[int, int, int]]) -> List[Tuple[int, int, int]]:
        """Merge overlapping tandem arrays, keeping the longest/strongest"""
        if not arrays:
            return []
        
        # Sort by start position
        sorted_arrays = sorted(arrays, key=lambda x: (x[0], -x[2]))  # Start pos, then by copy count desc
        merged = []
        
        for start, end, copies in sorted_arrays:
            # Check if this overlaps with any existing merged array
            overlaps = False
            for i, (m_start, m_end, m_copies) in enumerate(merged):
                if not (end <= m_start or start >= m_end):  # They overlap
                    # Keep the one with more copies, or if equal, the longer one
                    if copies > m_copies or (copies == m_copies and (end - start) > (m_end - m_start)):
                        merged[i] = (start, end, copies)
                    overlaps = True
                    break
            
            if not overlaps:
                merged.append((start, end, copies))
        
        return merged
    
    def score_cross_sequence_pattern(self, pattern: str, all_data: Dict[str, Dict[str, any]]) -> float:
        """
        Score a pattern based on cross-sequence presence and tandem strength
        """
        sequences_with_pattern = 0
        total_tandem_arrays = 0
        total_tandem_copies = 0
        max_single_array = 0
        total_occurrences = 0
        
        for seq_id, seq_data in all_data.items():
            if pattern in seq_data:
                sequences_with_pattern += 1
                data = seq_data[pattern]
                total_tandem_arrays += data['tandem_array_count']
                total_tandem_copies += data['total_tandem_copies']
                max_single_array = max(max_single_array, data['max_tandem_copies'])
                total_occurrences += data['total_occurrences']
        
        if sequences_with_pattern == 0:
            return 0
        
        # Scoring formula prioritizing:
        # 1. Cross-sequence presence
        # 2. Tandem array strength
        # 3. Pattern length
        # 4. Total occurrences
        
        sequence_coverage = sequences_with_pattern / len(self.sequences)
        pattern_length = len(pattern)
        avg_tandem_strength = total_tandem_copies / max(total_tandem_arrays, 1)
        
        score = (
            (sequence_coverage ** 1.2) *  # Strongly favor cross-sequence patterns
            (avg_tandem_strength ** 0.8) *  # Favor strong tandem arrays
            (pattern_length ** 0.2) *  # Modest preference for longer patterns
            math.log(total_occurrences + 1) *  # Logarithmic scaling for total occurrences
            math.log(max_single_array + 1)  # Bonus for very strong single arrays
        )
        
        return score

    def score_population_pattern(self, pattern: str, all_data: Dict) -> Tuple[float, float, float, float]:
        """
        Score a pattern for population mode (telomere read analysis).

        Scoring is based on per-read coverage fraction: how much of each positive
        read is explained by tandem copies of the pattern.  High coverage strongly
        indicates the read came from a t-circle.

        Returns (score, positive_read_fraction, mean_coverage, mean_max_copies)
        """
        pattern_length = len(pattern)
        total_reads = len(self.sequences)
        coverages: List[float] = []
        max_copies_list: List[int] = []

        for seq_id, seq_data in all_data.items():
            if pattern not in seq_data:
                continue
            data = seq_data[pattern]
            n_copies = data['non_overlapping_occurrences']
            if n_copies < 2:
                continue
            read_length = len(self.sequences[seq_id])
            coverage = min(1.0, n_copies * pattern_length / read_length)
            coverages.append(coverage)
            max_copies_list.append(data['max_tandem_copies'])

        if not coverages:
            return 0.0, 0.0, 0.0, 0.0

        positive_read_fraction = len(coverages) / total_reads
        mean_coverage = sum(coverages) / len(coverages)
        mean_max_copies = sum(max_copies_list) / len(max_copies_list)

        score = (
            (positive_read_fraction ** 0.8) *   # prevalence across reads
            (mean_coverage ** 1.5) *             # per-read coverage (key t-circle signal)
            ((mean_max_copies + 1) ** 0.4) *     # tandem strength
            (pattern_length ** 0.1)              # very mild length preference
        )
        return score, positive_read_fraction, mean_coverage, mean_max_copies

    def find_best_population_patterns(self, min_length: int, max_length: int,
                                      min_reads: int = 2, top_n: int = 15) -> List[Tuple[str, Dict, float]]:
        """
        Find best repeat patterns in population (read) mode.
        Scores patterns by per-read coverage rather than cross-sequence presence.
        Returns list of (pattern, stats, score) triples — same shape as
        find_best_cross_sequence_patterns so downstream code can be shared.
        Population-specific stats are stored in stats['population'].
        """
        self.pattern_file.write("="*90 + "\n")
        self.pattern_file.write(f"Pattern length range: {min_length}-{max_length} bp\n")
        self.pattern_file.write(f"Total reads: {len(self.sequences)}\n")
        self.pattern_file.write(f"Minimum reads for candidate: {min_reads}\n")

        all_data = self.analyze_all_sequences(min_length, max_length)

        all_patterns: Set[str] = set()
        for seq_data in all_data.values():
            all_patterns.update(seq_data.keys())
        self.pattern_file.write(f"\nFound {len(all_patterns)} unique patterns across all reads\n")

        scored_patterns: List[Tuple[str, Dict, float]] = []
        for pattern in all_patterns:
            reads_with_pattern = sum(
                1 for seq_data in all_data.values()
                if pattern in seq_data and seq_data[pattern]['non_overlapping_occurrences'] >= 2
            )
            if reads_with_pattern < min_reads:
                continue

            score, pos_frac, mean_cov, mean_copies = self.score_population_pattern(pattern, all_data)
            stats = self.compile_pattern_statistics(pattern, all_data)
            stats['population'] = {
                'positive_reads': reads_with_pattern,
                'total_reads': len(self.sequences),
                'positive_read_fraction': pos_frac,
                'mean_coverage': mean_cov,
                'mean_max_copies': mean_copies,
                'circle_confidence': pos_frac * mean_cov,
            }
            scored_patterns.append((pattern, stats, score))

        self.pattern_file.write(f"Found {len(scored_patterns)} candidate patterns in {min_reads}+ reads\n")

        # Reuse tandem-duplicate filter
        qualified = {p for p, _, _ in scored_patterns}
        scored_patterns = self._filter_tandem_duplicates(scored_patterns, qualified, min_length)
        self.pattern_file.write(f"After tandem-duplicate filtering: {len(scored_patterns)} patterns\n")

        scored_patterns.sort(key=lambda x: x[2], reverse=True)
        return scored_patterns[:top_n]

    def find_best_cross_sequence_patterns(self, min_length: int, max_length: int,
                                        min_sequences: int = 2, top_n: int = 20,
                                        scoring_weights: Dict[str, float] = None) -> List[Tuple[str, Dict, float]]:
        """
        Find the best patterns that appear across multiple sequences
        
        Args:
            scoring_weights: Dict to customize scoring behavior:
                - 'coverage_exp': Cross-sequence coverage weight (default 1.2)
                - 'tandem_exp': Tandem array strength weight (default 0.8)
                - 'length_exp': Pattern length weight (default 0.4)
                - 'use_log_length': Use logarithmic scaling for length (default False)
                
        Example scoring_weights for emphasizing longer patterns:
            {'length_exp': 0.8}  # Increase from 0.4 to 0.8
            {'length_exp': 1.0}  # Equal weight to cross-sequence coverage
            {'length_exp': 1.5}  # Strongly favor longer patterns
        """
        self.pattern_file.write("="*90+"\n")
        self.pattern_file.write(f"Pattern length range: {min_length}-{max_length} bp\n")
        self.pattern_file.write(f"Minimum sequences: {min_sequences}\n")
        
        if scoring_weights:
            self.pattern_file.write(f"Custom scoring weights: {scoring_weights}\n")
        
        # Analyze all sequences
        all_data = self.analyze_all_sequences(min_length, max_length)
        
        # Collect all unique patterns across sequences
        all_patterns = set()
        for seq_data in all_data.values():
            all_patterns.update(seq_data.keys())
        
        self.pattern_file.write(f"\nFound {len(all_patterns)} unique patterns across all sequences\n")
        
        # Score patterns and filter by minimum sequence requirement
        scored_patterns = []
        
        for pattern in all_patterns:
            sequences_with_pattern = sum(1 for seq_data in all_data.values() if pattern in seq_data)
            
            if sequences_with_pattern >= min_sequences:
            #    score = self.score_cross_sequence_pattern(pattern, all_data, scoring_weights)
                score = self.score_cross_sequence_pattern(pattern, all_data)
                
                # Compile comprehensive statistics
                stats = self.compile_pattern_statistics(pattern, all_data)
                scored_patterns.append((pattern, stats, score))
        
        self.pattern_file.write(f"Found {len(scored_patterns)} patterns in {min_sequences}+ sequences\n")

        # Filter out patterns that are exact tandem duplications of shorter primitives
        qualified_patterns = {p for p, _, _ in scored_patterns}
        scored_patterns = self._filter_tandem_duplicates(scored_patterns, qualified_patterns, min_length)
        self.pattern_file.write(f"After tandem-duplicate filtering: {len(scored_patterns)} patterns\n")

        # Sort by score and return top patterns
        scored_patterns.sort(key=lambda x: x[2], reverse=True)
        return scored_patterns[:top_n]
    
    def compile_pattern_statistics(self, pattern: str, all_data: Dict[str, Dict[str, any]]) -> Dict:
        """Compile comprehensive statistics for a pattern"""
        stats = {
            'pattern': pattern,
            'pattern_length': len(pattern),
            'sequences_present': [],
            'sequence_coverage': 0,
            'total_occurrences': 0,
            'total_tandem_arrays': 0,
            'total_tandem_copies': 0,
            'max_tandem_array': 0,
            'avg_tandem_strength': 0,
            'sequence_details': {}
        }
        
        for seq_id, seq_data in all_data.items():
            if pattern in seq_data:
                data = seq_data[pattern]
                stats['sequences_present'].append(seq_id)
                stats['total_occurrences'] += data['total_occurrences']
                stats['total_tandem_arrays'] += data['tandem_array_count']
                stats['total_tandem_copies'] += data['total_tandem_copies']
                stats['max_tandem_array'] = max(stats['max_tandem_array'], data['max_tandem_copies'])
                
                stats['sequence_details'][seq_id] = {
                    'occurrences': data['total_occurrences'],
                    'tandem_arrays': data['tandem_array_count'],
                    'tandem_copies': data['total_tandem_copies'],
                    'max_array_size': data['max_tandem_copies'],
                    'array_positions': data['tandem_arrays']
                }
        
        stats['sequence_coverage'] = len(stats['sequences_present']) / len(self.sequences)
        stats['avg_tandem_strength'] = stats['total_tandem_copies'] / max(stats['total_tandem_arrays'], 1)
        
        return stats
    
    def visualize_pattern_across_sequences(self, pattern: str, stats: Dict, context: int = 50):
        """
        Visualize how a pattern appears across sequences
        """
        self.pattern_file.write(f"\n{'='*90}\n")
        self.pattern_file.write(f"PATTERN VISUALIZATION: {pattern[:40]}{'...' if len(pattern) > 40 else ''}\n")
        self.pattern_file.write(f"{'='*90}\n")
        self.pattern_file.write(f"Length: {stats['pattern_length']} bp\n")
        self.pattern_file.write(f"Present in: {len(stats['sequences_present'])}/{len(self.sequences)} sequences\n")
        self.pattern_file.write(f"Total tandem arrays: {stats['total_tandem_arrays']}\n")
        self.pattern_file.write(f"Total tandem copies: {stats['total_tandem_copies']}\n")
        self.pattern_file.write(f"Largest array: {stats['max_tandem_array']} copies\n")
        
        for seq_id in stats['sequences_present']:
            self.pattern_file.write(f"\n--- {seq_id} ---\n")
            seq_details = stats['sequence_details'][seq_id]
            self.pattern_file.write(f"Arrays: {seq_details['tandem_arrays']}, Max size: {seq_details['max_array_size']} copies\n")
            
            sequence = self.sequences[seq_id]
            
            for start, end, copies in seq_details['array_positions']:
                array_length = end - start
                self.pattern_file.write(f"  Array: pos {start:,}-{end:,} ({array_length:,} bp, {copies} copies)\n")
                
                # Show context
                context_start = max(0, start - context)
                context_end = min(len(sequence), end + context)
                
                before = sequence[context_start:start]
                array_seq = sequence[start:end]
                after = sequence[end:context_end]
                
                self.pattern_file.write(f"    ...{before[-30:]}[{array_seq[:60]}{'...' if len(array_seq) > 60 else ''}]{after[:30]}...\n")


def analyze_chromosome_ends(sequences: Dict[str, str], pattern_file: TextIOWrapper, min_pattern: int, max_pattern: int,
                          scoring_weights: Dict[str, float] = None):
    """
    Main analysis function for chromosome end repeat patterns
    
    Args:
        sequences: Dict of sequence_id -> sequence_string
        min_pattern: Minimum pattern length
        max_pattern: Maximum pattern length
        scoring_weights: Optional dict to customize scoring:
            - 'coverage_exp': Weight for cross-sequence coverage (default 1.2)
            - 'tandem_exp': Weight for tandem strength (default 0.8)
            - 'length_exp': Weight for pattern length (default 0.4)
            - 'use_log_length': Use log scaling for length (default False)
            
    Example usage:
        # Default: modest preference for longer patterns
        analyze_chromosome_ends(seqs)
        
        # Strong preference for longer patterns
        analyze_chromosome_ends(seqs, scoring_weights={'length_exp': 1.2})
        
        # Extreme preference for longer patterns
        analyze_chromosome_ends(seqs, scoring_weights={'length_exp': 2.0})
        
        # Balanced: equal weight to all factors
        analyze_chromosome_ends(seqs, scoring_weights={
            'coverage_exp': 1.0,
            'tandem_exp': 1.0,
            'length_exp': 1.0
        })
    """
    pattern_file.write("CHROMOSOME END REPEAT ANALYSIS\n")
    pattern_file.write("="*90 + "\n")
    pattern_file.write(f"  Pattern length range : {min_pattern}–{max_pattern} bp\n")
    pattern_file.write(f"  Sequences            : {len(sequences)}\n")
    pattern_file.write("="*90 + "\n\n")
    pattern_file.write(f"{'Sequence':<8}  {'Length':>10}\n")
    pattern_file.write("-"*25 + "\n")

    finder = ChromosomeEndRepeatFinder(sequences, pattern_file)
    pattern_file.write("\n" + "="*90 + "\n")
    pattern_file.write("PER-SEQUENCE ANALYSIS\n")
    pattern_file.write("="*90 + "\n")

    # Find best cross-sequence patterns
    results = finder.find_best_cross_sequence_patterns(
        min_length=min_pattern,
        max_length=max_pattern,
        min_sequences=2,
        top_n=15,
        scoring_weights=scoring_weights
    )
    
    if not results:
        pattern_file.write("No cross-sequence patterns found!\n")
        return finder, None, results

    n_seqs = len(finder.sequences)

    # ── Top candidates table ───────────────────────────────────────────────
    pattern_file.write(f"\n{'='*90}\n")
    pattern_file.write("TOP CANDIDATE PATTERNS\n")
    pattern_file.write("="*90 + "\n")
    pattern_file.write(f"{'Rank':>4}  {'Length':>6}  {'Seqs':>5}  {'Arrays':>6}  {'Copies':>6}  {'MaxArr':>6}  {'Score':>8}  Pattern Preview\n")
    pattern_file.write("-"*90 + "\n")

    for i, (pattern, stats, score) in enumerate(results):
        preview = pattern[:40] + "..." if len(pattern) > 40 else pattern
        seqs = len(stats['sequences_present'])
        pattern_file.write(
            f"{i+1:4d}  {stats['pattern_length']:6d}  {seqs:4d}/{n_seqs:<3d}"
            f"  {stats['total_tandem_arrays']:6d}  {stats['total_tandem_copies']:6d}"
            f"  {stats['max_tandem_array']:6d}  {score:8.2f}  {preview}\n"
        )

    # ── Best pattern summary ───────────────────────────────────────────────
    best_pattern, best_stats, best_score = results[0]
    best_seqs = len(best_stats['sequences_present'])

    pattern_file.write(f"\n{'='*90}\n")
    pattern_file.write("BEST PATTERN\n")
    pattern_file.write("="*90 + "\n")
    pattern_file.write(f"  Length  : {best_stats['pattern_length']} bp\n")
    pattern_file.write(f"  Present : {best_seqs}/{n_seqs} sequences\n")
    pattern_file.write(f"  Arrays  : {best_stats['total_tandem_arrays']} total, largest {best_stats['max_tandem_array']} copies\n")
    pattern_file.write(f"  Score   : {best_score:.2f}\n")
    pattern_file.write(f"  Sequence:\n    {best_pattern}\n")

    # Print summary to screen
    print(f"\n{'='*70}")
    print(f"BEST PATTERN FOUND")
    print(f"{'='*70}")
    print(f"  Length  : {best_stats['pattern_length']} bp")
    print(f"  Present : {best_seqs}/{n_seqs} sequences")
    print(f"  Arrays  : {best_stats['total_tandem_arrays']} total, largest {best_stats['max_tandem_array']} copies")
    print(f"  Score   : {best_score:.2f}")
    print(f"  Sequence:")
    print(f"    {best_pattern}")
    print(f"{'='*70}\n")

    # ── Per-sequence visualization for best pattern ────────────────────────
    finder.visualize_pattern_across_sequences(best_pattern, best_stats)

    return finder, best_pattern, results


# Create test data that matches your requirements
def create_test_chromosome_sequences():
    """Create test sequences with large patterns (100-200bp) that form tandem repeats"""
    
    # Large repeat pattern (~150bp) - like a telomeric repeat unit
    large_pattern = (
        "TTAGGGTTACGTACGTACGTACGTACGTACGTACGTACGT"  # 40bp
        "AAAAAGGGGCCCCTTTTTATCGATCGATCGATCGATCGAT"  # 40bp  
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"  # 40bp
        "NNNATCGATCG"  # 11bp - total = 131bp
    )
    
    # Alternative pattern (~120bp)
    alt_pattern = (
        "CCCCGGGGAAAAATTTTTATCGATCGATCGATCGATCGAT"  # 40bp
        "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"  # 40bp
        "TTTAAACCCGGGATCGATCGATCGATCG"  # 27bp - total = 107bp
    )
    
    sequences = {}
    
    # Chromosome 1p - strong tandem repeat of large pattern
    seq1 = "N" * 200
    seq1 += large_pattern * 5  # 5 consecutive copies (tandem array)
    seq1 += "ATCGATCGATCG" * 20  # spacer
    seq1 += large_pattern * 3  # 3 more consecutive copies
    seq1 += "N" * 150
    sequences["chr1p_telomere"] = seq1
    
    # Chromosome 2q - same pattern, different arrangement  
    seq2 = "N" * 180
    seq2 += large_pattern * 7  # 7 consecutive copies
    seq2 += "GCTAGCTA" * 25  # different spacer
    seq2 += large_pattern * 2  # 2 more copies
    seq2 += "N" * 100
    sequences["chr2q_telomere"] = seq2
    
    # Chromosome 3p - rotated version of the pattern
    rotated_pattern = large_pattern[50:] + large_pattern[:50]  # Rotate by 50bp
    seq3 = "N" * 220
    seq3 += rotated_pattern * 4  # 4 consecutive copies of rotated pattern
    seq3 += "AAAATTTT" * 30  # spacer
    seq3 += rotated_pattern * 6  # 6 more copies
    seq3 += "N" * 80
    sequences["chr3p_telomere"] = seq3
    
    # Chromosome 4q - alternative pattern
    seq4 = "N" * 160
    seq4 += alt_pattern * 3  # 3 copies of different pattern
    seq4 += "CCCCGGGG" * 15  # spacer
    seq4 += alt_pattern * 4  # 4 more copies
    seq4 += "N" * 200
    sequences["chr4q_telomere"] = seq4
    
    # Chromosome 5p - both patterns
    seq5 = "N" * 300
    seq5 += large_pattern * 3  # 3 copies of main pattern
    seq5 += "TTTTAAAA" * 25  # spacer
    seq5 += alt_pattern * 2   # 2 copies of alt pattern
    seq5 += "GGGGCCCC" * 10   # spacer
    seq5 += large_pattern * 4  # 4 more copies of main pattern
    seq5 += "N" * 120
    sequences["chr5p_telomere"] = seq5
    
    return sequences, large_pattern, alt_pattern

def _create_stats_file(output_file) -> TextIOWrapper:
    base, ext = splitext(output_file)
    pattern_file_name = f"{base}_pattern.txt"
    pattern_output_file = open(pattern_file_name, 'w')
    return pattern_output_file

def analyze_population_reads(sequences: Dict[str, str], pattern_file: TextIOWrapper,
                             min_pattern: int, max_pattern: int,
                             total_raw_reads: Optional[int] = None,
                             max_workers: Optional[int] = None):
    """
    Top-level function for population mode: find repeat patterns and assess
    whether a t-circle is present in the read set.

    Analogous to analyze_chromosome_ends() but uses per-read coverage scoring
    instead of cross-sequence presence scoring, and outputs a circle confidence score.
    """
    pattern_file.write("POPULATION TELOMERE READ ANALYSIS\n")
    pattern_file.write("="*90 + "\n")
    pattern_file.write(f"  Pattern length range : {min_pattern}–{max_pattern} bp\n")
    if total_raw_reads is not None:
        pattern_file.write(f"  Total raw reads      : {total_raw_reads}\n")
        pattern_file.write(f"  Telomeric reads      : {len(sequences)}\n")
    else:
        pattern_file.write(f"  Total reads          : {len(sequences)}\n")
    pattern_file.write("="*90 + "\n\n")

    finder = ChromosomeEndRepeatFinder(sequences, pattern_file, max_workers=max_workers)
    pattern_file.write("\n" + "="*90 + "\n")
    pattern_file.write("PER-READ ANALYSIS\n")
    pattern_file.write("="*90 + "\n")

    results = finder.find_best_population_patterns(
        min_length=min_pattern,
        max_length=max_pattern,
        min_reads=2,
        top_n=15,
    )

    if not results:
        pattern_file.write("No repeated patterns found across reads.\n")
        print("No repeated patterns found in the reads.")
        return finder, None, results

    n_reads = len(finder.sequences)

    # ── Top candidates table ───────────────────────────────────────────────
    pattern_file.write(f"\n{'='*90}\n")
    pattern_file.write("TOP CANDIDATE PATTERNS\n")
    pattern_file.write("="*90 + "\n")
    pattern_file.write(
        f"{'Rank':>4}  {'Length':>6}  {'Reads':>9}  "
        f"{'MeanCov':>8}  {'MeanCop':>8}  {'CircConf':>9}  {'Score':>8}  Pattern Preview\n"
    )
    pattern_file.write("-"*90 + "\n")

    for i, (pattern, stats, score) in enumerate(results):
        pop = stats['population']
        preview = pattern[:40] + "..." if len(pattern) > 40 else pattern
        pattern_file.write(
            f"{i+1:4d}  {stats['pattern_length']:6d}  "
            f"{pop['positive_reads']:4d}/{n_reads:<4d}  "
            f"{pop['mean_coverage']:8.1%}  "
            f"{pop['mean_max_copies']:8.1f}  "
            f"{pop['circle_confidence']:9.3f}  "
            f"{score:8.2f}  {preview}\n"
        )

    # ── Circle detection verdict ───────────────────────────────────────────
    best_pattern, best_stats, best_score = results[0]
    pop = best_stats['population']
    conf = pop['circle_confidence']

    if conf >= 0.4:
        verdict = "STRONG CIRCLE SIGNAL"
    elif conf >= 0.2:
        verdict = "MODERATE CIRCLE SIGNAL"
    elif conf >= 0.1:
        verdict = "WEAK / AMBIGUOUS SIGNAL"
    else:
        verdict = "NO CIRCLE DETECTED"

    pattern_file.write(f"\n{'='*90}\n")
    pattern_file.write("CIRCLE DETECTION RESULT\n")
    pattern_file.write("="*90 + "\n")
    pattern_file.write(f"  Verdict            : {verdict}\n")
    pattern_file.write(f"  Circle confidence  : {conf:.3f}\n")
    pattern_file.write(f"  Reads with pattern : {pop['positive_reads']}/{n_reads} ({pop['positive_read_fraction']:.1%})\n")
    pattern_file.write(f"  Mean read coverage : {pop['mean_coverage']:.1%}\n")
    pattern_file.write(f"  Mean tandem copies : {pop['mean_max_copies']:.1f}\n")
    pattern_file.write(f"\n{'='*90}\n")
    pattern_file.write("BEST PATTERN\n")
    pattern_file.write("="*90 + "\n")
    pattern_file.write(f"  Length  : {best_stats['pattern_length']} bp\n")
    pattern_file.write(f"  Score   : {best_score:.2f}\n")
    pattern_file.write(f"  Sequence:\n    {best_pattern}\n")

    # Print to screen
    print(f"\n{'='*70}")
    print(f"CIRCLE DETECTION: {verdict}")
    print(f"{'='*70}")
    print(f"  Confidence  : {conf:.3f}")
    print(f"  Reads       : {pop['positive_reads']}/{n_reads} ({pop['positive_read_fraction']:.1%})")
    print(f"  Coverage    : {pop['mean_coverage']:.1%} mean per positive read")
    print(f"  Copies      : {pop['mean_max_copies']:.1f} mean tandem copies per positive read")
    print(f"  Pattern     : {best_stats['pattern_length']} bp")
    print(f"    {best_pattern}")
    print(f"{'='*70}\n")

    finder.visualize_pattern_across_sequences(best_pattern, best_stats)
    return finder, best_pattern, results


def pattern_finder_execute(telomeres: List[Optional[TelomereSequence]], config: Config) -> Optional[str]:
    pattern_file = _create_stats_file(config.output_file)
    test_dict: Dict[str, str] = {}
    for telomer in telomeres:
        if telomer and telomer.sequence:
            test_dict[telomer.chromosome_end_id] = telomer.sequence

    finder, best_pattern, results = analyze_chromosome_ends(
        test_dict, pattern_file, config.min_pattern_length, config.max_pattern_length
    )

    if best_pattern:
        return best_pattern
    return None


def pattern_finder_execute_population(config: Config) -> Optional[str]:
    """
    Entry point for population mode pattern finding.
    Extracts telomeric regions from both ends of each read, converts all to
    TG orientation, then runs population scoring and prints a circle detection verdict.
    """
    from data_io.fasta_reader import FastaReader

    pattern_file = _create_stats_file(config.output_file)

    max_workers = _resolve_workers(config.workers)
    if max_workers:
        print(f"  Using {max_workers} worker processes")

    reader = FastaReader(config.fasta_file_path, max_ends=0)
    sequences, total_raw = reader.parse_fasta_population(
        threshold=config.telomere_threshold, max_workers=max_workers
    )

    if not sequences:
        raise ValueError("No telomeric sequences found in FASTA file")

    # Write extracted telomeric reads to a FASTA file
    base, ext = splitext(config.output_file)
    telomere_fasta_path = f"{base}_telomeres.fasta"
    with open(telomere_fasta_path, 'w') as tf:
        for read_id, seq in sequences.items():
            tf.write(f">{read_id}\n{seq}\n")
    print(f"  Telomere extraction: {len(sequences)} telomeric reads from {total_raw} total reads")
    print(f"  Telomeric reads written to: {telomere_fasta_path}")

    _, best_pattern, _ = analyze_population_reads(
        sequences, pattern_file, config.min_pattern_length, config.max_pattern_length,
        total_raw_reads=total_raw, max_workers=max_workers,
    )
    return best_pattern



# Test the system
if __name__ == "__main__":
    print("Creating test chromosome end sequences...")
    test_sequences, expected_main, expected_alt = create_test_chromosome_sequences()
    
    print(f"\nExpected main pattern ({len(expected_main)} bp):")
    print(f"  {expected_main[:60]}...")
    
    print(f"\nExpected alt pattern ({len(expected_alt)} bp):")  
    print(f"  {expected_alt[:60]}...")
    
    # Run analysis
    finder, best_pattern, results = analyze_chromosome_ends(test_sequences, 100, 180)
    
    if results:
        print(f"\n{'='*90}")
        print("VALIDATION")
        print("="*90)
        
        best_pattern_found = results[0][0]
        print(f"Best pattern found ({len(best_pattern_found)} bp):")
        print(f"  {best_pattern_found[:60]}...")
        
        # Check if found pattern matches expected (considering rotations)
        finder_obj = ChromosomeEndRepeatFinder(test_sequences)
        expected_canonical = finder_obj.get_canonical_rotation(expected_main)
        found_canonical = finder_obj.get_canonical_rotation(best_pattern_found)
        
        print(f"\nCanonical forms match: {expected_canonical == found_canonical}")
        if expected_canonical != found_canonical:
            alt_canonical = finder_obj.get_canonical_rotation(expected_alt)
            print(f"Matches alt pattern: {alt_canonical == found_canonical}")