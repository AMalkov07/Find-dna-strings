import math
from collections import defaultdict, Counter
from typing import List, Tuple, Dict, Set, Optional
from utils.data_structures import TelomereSequence

class ChromosomeEndRepeatFinder:
    def __init__(self, sequences: Dict[str, str]):
        """
        Initialize with multiple DNA sequences (e.g., chromosome ends)
        sequences: Dict mapping sequence_id -> sequence_string
        """
        self.sequences = {seq_id: seq.upper().replace('N', '') for seq_id, seq in sequences.items()}
        self.sequence_ids = list(self.sequences.keys())
        print(f"Loaded {len(self.sequences)} sequences")
        for seq_id, seq in self.sequences.items():
            print(f"  {seq_id}: {len(seq):,} bp")
    
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
        
        rotations = []
        for i in range(len(pattern)):
            rotation = pattern[i:] + pattern[:i]
            rotations.append(rotation)
        
        return min(rotations)
    
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
        print(f"    Consolidating rotational patterns from {len(raw_patterns)} raw patterns...")
        
        # Group patterns by their canonical rotation
        canonical_groups = defaultdict(list)  # canonical -> [(pattern, positions), ...]
        
        for pattern, position_set in raw_patterns.items():
            canonical = self.get_canonical_rotation(pattern)
            canonical_groups[canonical].append((pattern, position_set))
        
        print(f"    Grouped into {len(canonical_groups)} canonical pattern families")
        
        # For each canonical group, find the best representative and all its positions
        final_patterns = {}
        
        for canonical, pattern_list in canonical_groups.items():
            # Find which rotation has the most occurrences
            best_pattern = canonical
            best_positions = set()
            max_occurrences = 0
            
            for pattern, positions in pattern_list:
                if len(positions) > max_occurrences:
                    max_occurrences = len(positions)
                    best_pattern = pattern
                    best_positions = positions.copy()
                elif len(positions) == max_occurrences:
                    # Same number of occurrences, combine positions
                    best_positions.update(positions)
            
            # Now find ALL occurrences of this pattern and its rotations efficiently
            all_rotation_positions = self.find_all_rotations_efficiently(sequence, best_pattern)
            
            # Combine with known positions
            all_rotation_positions.update(best_positions)
            
            final_patterns[canonical] = sorted(list(all_rotation_positions))
        
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
        print(f"  Building suffix array for {len(sequence):,} bp sequence...")
        suffix_array = self._build_suffix_array(sequence)
        
        print(f"  Building LCP array...")
        lcp_array = self._build_lcp_array(sequence, suffix_array)
        
        print(f"  Extracting repeated patterns {min_length}-{max_length} bp...")
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
        
        print(f"    Found {len(raw_patterns)} raw patterns before rotation consolidation")
        
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
        Analyze all sequences to find repeated patterns and their tandem arrays
        SUPER OPTIMIZED VERSION - handles rotations efficiently upfront
        Returns: Dict[seq_id -> Dict[pattern -> analysis_data]]
        """
        all_sequence_data = {}
        
        for seq_id, sequence in self.sequences.items():
            print(f"\nAnalyzing sequence {seq_id}...")
            
            # Step 1: Find all repeated patterns with rotation consolidation built-in
            # This now returns canonical patterns with ALL positions (including rotations)
            repeated_patterns = self.find_repeated_patterns_in_sequence(sequence, min_length, max_length)
            print(f"  Found {len(repeated_patterns)} consolidated patterns")
            
            # Step 2: For each consolidated pattern, find tandem arrays using all positions
            sequence_analysis = {}
            
            for canonical_pattern, all_positions in repeated_patterns.items():
                if len(all_positions) >= 2:  # Must appear at least twice
                    pattern_len = len(canonical_pattern)
                    
                    # Get non-overlapping positions for proper counting
                    non_overlapping_positions = self.get_non_overlapping_positions(all_positions, pattern_len)
                    
                    print(f"    Processing pattern of length {pattern_len}: "
                          f"{len(all_positions)} total / {len(non_overlapping_positions)} non-overlapping occurrences")
                    
                    # Find tandem arrays using all known positions (including overlapping)
                    # The tandem detection itself handles overlap correctly
                    tandem_arrays = self.find_tandem_arrays_optimized(all_positions, pattern_len)
                    
                    sequence_analysis[canonical_pattern] = {
                        'positions': all_positions,  # All positions including rotations and overlaps
                        'non_overlapping_positions': non_overlapping_positions,  # Non-overlapping only
                        'total_occurrences': len(all_positions),  # Total including overlaps
                        'non_overlapping_occurrences': len(non_overlapping_positions),  # Non-overlapping count
                        'tandem_arrays': tandem_arrays,
                        'tandem_array_count': len(tandem_arrays),
                        'total_tandem_copies': sum(copies for _, _, copies in tandem_arrays),
                        'max_tandem_copies': max((copies for _, _, copies in tandem_arrays), default=0),
                        'pattern_length': len(canonical_pattern)
                    }
            
            all_sequence_data[seq_id] = sequence_analysis
            print(f"  Found {len(sequence_analysis)} patterns with tandem arrays")
        
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
        
        #score = (
            #(sequence_coverage ** 1.2) *  # Strongly favor cross-sequence patterns
            #(avg_tandem_strength ** 0.4) *  # Favor strong tandem arrays
            #(pattern_length ** 100) *  # Modest preference for longer patterns
            #math.log(total_occurrences + 1) *  # Logarithmic scaling for total occurrences
            #math.log(max_single_array + 1)  # Bonus for very strong single arrays
        #)
        score = (
            (pattern_length)  # Modest preference for longer patterns
        )
        
        return score
    
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
        print("="*90)
        print("CROSS-SEQUENCE REPEAT PATTERN ANALYSIS")
        print("="*90)
        print(f"Pattern length range: {min_length}-{max_length} bp")
        print(f"Minimum sequences: {min_sequences}")
        
        if scoring_weights:
            print(f"Custom scoring weights: {scoring_weights}")
        
        # Analyze all sequences
        all_data = self.analyze_all_sequences(min_length, max_length)
        
        # Collect all unique patterns across sequences
        all_patterns = set()
        for seq_data in all_data.values():
            all_patterns.update(seq_data.keys())
        
        print(f"\nFound {len(all_patterns)} unique patterns across all sequences")
        
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
        
        print(f"Found {len(scored_patterns)} patterns in {min_sequences}+ sequences")
        
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
        print(f"\n{'='*90}")
        print(f"PATTERN VISUALIZATION: {pattern[:40]}{'...' if len(pattern) > 40 else ''}")
        print(f"{'='*90}")
        print(f"Length: {stats['pattern_length']} bp")
        print(f"Present in: {len(stats['sequences_present'])}/{len(self.sequences)} sequences")
        print(f"Total tandem arrays: {stats['total_tandem_arrays']}")
        print(f"Total tandem copies: {stats['total_tandem_copies']}")
        print(f"Largest array: {stats['max_tandem_array']} copies")
        
        for seq_id in stats['sequences_present']:
            print(f"\n--- {seq_id} ---")
            seq_details = stats['sequence_details'][seq_id]
            print(f"Arrays: {seq_details['tandem_arrays']}, Max size: {seq_details['max_array_size']} copies")
            
            sequence = self.sequences[seq_id]
            
            for start, end, copies in seq_details['array_positions']:
                array_length = end - start
                print(f"  Array: pos {start:,}-{end:,} ({array_length:,} bp, {copies} copies)")
                
                # Show context
                context_start = max(0, start - context)
                context_end = min(len(sequence), end + context)
                
                before = sequence[context_start:start]
                array_seq = sequence[start:end]
                after = sequence[end:context_end]
                
                print(f"    ...{before[-30:]}[{array_seq[:60]}{'...' if len(array_seq) > 60 else ''}]{after[:30]}...")


def analyze_chromosome_ends(sequences: Dict[str, str], min_pattern: int, max_pattern: int,
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
    print("CHROMOSOME END REPEAT ANALYSIS")
    print("="*90)
    
    finder = ChromosomeEndRepeatFinder(sequences)
    
    # Find best cross-sequence patterns
    results = finder.find_best_cross_sequence_patterns(
        min_length=min_pattern,
        max_length=max_pattern,
        min_sequences=2,
        top_n=15,
        scoring_weights=scoring_weights
    )
    
    if not results:
        print("No cross-sequence patterns found!")
        return finder, None, results
    
    # Display results table
    print(f"\n{'='*90}")
    print("TOP CROSS-SEQUENCE PATTERNS")
    print("="*90)
    print(f"{'Rank':>4} {'Length':>6} {'Seqs':>4} {'Arrays':>6} {'Copies':>7} {'Max Array':>9} {'Score':>8} {'Pattern Preview':>30}")
    print("-"*90)
    
    for i, (pattern, stats, score) in enumerate(results):
        preview = pattern[:25] + "..." if len(pattern) > 25 else pattern
        print(f"{i+1:4d} {stats['pattern_length']:6d} {len(stats['sequences_present']):4d} "
              f"{stats['total_tandem_arrays']:6d} {stats['total_tandem_copies']:7d} "
              f"{stats['max_tandem_array']:9d} {score:8.2f} {preview:>30}")
    
    # Detailed analysis of top pattern
    best_pattern, best_stats, best_score = results[0]
    finder.visualize_pattern_across_sequences(best_pattern, best_stats)

    print(f"best_pattern: {best_pattern}")
    for r in results:
        print(r)
    
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


def execute(telomeres: List[Optional[TelomereSequence]]) -> Optional[str]:
    test_dict: Dict[str, str] = {}
    for telomer in telomeres:
        if telomer and telomer.sequence:
            test_dict[telomer.chromosome_end_id] = telomer.sequence

    finder, best_pattern, results = analyze_chromosome_ends(test_dict, 50, 300)

    if best_pattern:
      return best_pattern
    return None



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