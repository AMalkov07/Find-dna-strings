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
    
    def get_all_rotations(self, pattern: str) -> List[str]:
        """Get all circular rotations of a pattern"""
        rotations = []
        for i in range(len(pattern)):
            rotation = pattern[i:] + pattern[:i]
            rotations.append(rotation)
        return rotations
    
    def is_rotation(self, str1: str, str2: str) -> bool:
        """
        Fast check if str2 is a rotation of str1 using the doubling trick
        """
        if len(str1) != len(str2):
            return False
        if str1 == str2:
            return True
        return str2 in (str1 + str1)
    
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
    
    def find_all_non_overlapping_occurrences(self, sequence: str, pattern: str) -> List[int]:
        """
        Find all non-overlapping occurrences of a pattern (including rotations) in a sequence
        This is the key function that ensures we only work with non-overlapping positions
        """
        pattern_len = len(pattern)
        doubled_pattern = pattern + pattern
        
        # Find all positions (including overlapping) first
        all_positions = []
        for i in range(len(sequence) - pattern_len + 1):
            candidate = sequence[i:i + pattern_len]
            if candidate in doubled_pattern:  # Checks all rotations via doubling trick
                all_positions.append(i)
        
        # Filter to non-overlapping
        return self.get_non_overlapping_positions(all_positions, pattern_len)
    
    def find_repeated_patterns_in_sequence(self, sequence: str, min_length: int = 100, 
                                         max_length: int = 200, step_size: int = 1) -> Dict[str, List[int]]:
        """
        Use suffix array + LCP to find repeated patterns, returning NON-OVERLAPPING positions only
        Returns: Dict[canonical_pattern -> List[non_overlapping_positions]]
        """
        print(f"  Building suffix array for {len(sequence):,} bp sequence...")
        suffix_array = self._build_suffix_array(sequence)
        
        print(f"  Building LCP array...")
        lcp_array = self._build_lcp_array(sequence, suffix_array)
        
        print(f"  Extracting repeated patterns {min_length}-{max_length} bp (step={step_size})...")
        
        # Track LCP statistics
        max_lcp_found = 0
        lcp_distribution = defaultdict(int)
        
        # First pass: collect all patterns from LCP (may include overlapping positions)
        raw_patterns = defaultdict(set)
        
        for i, lcp_len in enumerate(lcp_array):
            if lcp_len > max_lcp_found:
                max_lcp_found = lcp_len
            
            if lcp_len >= min_length:
                lcp_distribution[min(lcp_len // 10 * 10, max_length)] += 1
                
                pos1 = suffix_array[i]
                pos2 = suffix_array[i + 1]
                
                max_extract = min(lcp_len, max_length)
                
                for length in range(min_length, max_extract + 1, step_size):
                    if pos1 + length <= len(sequence):
                        pattern = sequence[pos1:pos1 + length]
                        canonical = self.get_canonical_rotation(pattern)
                        raw_patterns[canonical].add(pos1)
                        raw_patterns[canonical].add(pos2)
        
        print(f"    Max LCP found: {max_lcp_found} bp")
        print(f"    LCP length distribution:")
        for length_bucket in sorted(lcp_distribution.keys()):
            print(f"      {length_bucket}-{length_bucket+9}bp: {lcp_distribution[length_bucket]} patterns")
        print(f"    Found {len(raw_patterns)} unique canonical patterns from LCP")
        
        # Second pass: for each canonical pattern, find ALL non-overlapping occurrences
        # (including all rotations) using comprehensive search
        final_patterns = {}
        
        print(f"  Finding all non-overlapping occurrences for each pattern...")
        for canonical_pattern, lcp_positions in raw_patterns.items():
            if len(lcp_positions) >= 2:  # Must appear at least twice
                # Find ALL non-overlapping occurrences including all rotations
                non_overlapping_positions = self.find_all_non_overlapping_occurrences(
                    sequence, canonical_pattern
                )
                
                if len(non_overlapping_positions) >= 2:
                    final_patterns[canonical_pattern] = non_overlapping_positions
        
        print(f"  Final: {len(final_patterns)} patterns with 2+ non-overlapping occurrences")
        
        # Show pattern length distribution
        length_dist = defaultdict(int)
        for pattern in final_patterns.keys():
            length_bucket = len(pattern) // 10 * 10
            length_dist[length_bucket] += 1
        
        if length_dist:
            print(f"  Pattern length distribution:")
            for length_bucket in sorted(length_dist.keys()):
                print(f"    {length_bucket}-{length_bucket+9}bp: {length_dist[length_bucket]} patterns")
        
        return final_patterns
    
    def find_tandem_arrays_optimized(self, known_positions: List[int], pattern_len: int) -> List[Tuple[int, int, int]]:
        """
        Find tandem arrays from known non-overlapping positions
        Args:
            known_positions: List of NON-OVERLAPPING positions where pattern occurs
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
    
    def analyze_all_sequences(self, min_length: int = 100, max_length: int = 200, 
                            step_size: int = 1) -> Dict[str, Dict[str, any]]:
        """
        Analyze all sequences to find repeated patterns and their tandem arrays
        All occurrence counts are NON-OVERLAPPING throughout
        """
        all_sequence_data = {}
        
        for seq_id, sequence in self.sequences.items():
            print(f"\nAnalyzing sequence {seq_id}...")
            
            # Find all repeated patterns with NON-OVERLAPPING positions only
            repeated_patterns = self.find_repeated_patterns_in_sequence(
                sequence, min_length, max_length, step_size
            )
            print(f"  Found {len(repeated_patterns)} patterns")
            
            # For each pattern, find tandem arrays
            sequence_analysis = {}
            
            for canonical_pattern, non_overlapping_positions in repeated_patterns.items():
                pattern_len = len(canonical_pattern)
                
                print(f"    Pattern length {pattern_len}: {len(non_overlapping_positions)} non-overlapping occurrences")
                
                # Find tandem arrays using non-overlapping positions
                tandem_arrays = self.find_tandem_arrays_optimized(non_overlapping_positions, pattern_len)
                
                sequence_analysis[canonical_pattern] = {
                    'positions': non_overlapping_positions,  # NON-OVERLAPPING positions
                    'occurrences': len(non_overlapping_positions),  # Count of non-overlapping occurrences
                    'tandem_arrays': tandem_arrays,
                    'tandem_array_count': len(tandem_arrays),
                    'total_tandem_copies': sum(copies for _, _, copies in tandem_arrays),
                    'max_tandem_copies': max((copies for _, _, copies in tandem_arrays), default=0),
                    'pattern_length': len(canonical_pattern)
                }
            
            all_sequence_data[seq_id] = sequence_analysis
            print(f"  Analyzed {len(sequence_analysis)} patterns")
        
        return all_sequence_data
    
    def score_cross_sequence_pattern(self, pattern, all_data, scoring_weights=None):
        """
        Score a pattern based on cross-sequence presence and tandem strength
        All occurrence counts are non-overlapping
        """
        # Default weights
        if scoring_weights is None:
            scoring_weights = {
                'coverage_exp': 1.2,
                'tandem_exp': 0.8,
                'length_exp': 0.4,
                'use_log_length': False
            }
        
        sequences_with_pattern = 0
        total_tandem_arrays = 0
        total_tandem_copies = 0
        max_single_array = 0
        total_occurrences = 0  # Non-overlapping occurrences
        
        canonical_pattern = self.get_canonical_rotation(pattern)
        
        for seq_id, seq_data in all_data.items():
            if canonical_pattern in seq_data:
                sequences_with_pattern += 1
                data = seq_data[canonical_pattern]
                total_tandem_arrays += data['tandem_array_count']
                total_tandem_copies += data['total_tandem_copies']
                max_single_array = max(max_single_array, data['max_tandem_copies'])
                total_occurrences += data['occurrences']  # Already non-overlapping
        
        if sequences_with_pattern == 0:
            return 0
        
        # Calculate components
        sequence_coverage = sequences_with_pattern / len(self.sequences)
        pattern_length = len(pattern)
        avg_tandem_strength = total_tandem_copies / max(total_tandem_arrays, 1)
        
        # Apply length scaling
        if scoring_weights.get('use_log_length', False):
            length_component = math.log(pattern_length + 1) ** scoring_weights['length_exp']
        else:
            length_component = pattern_length ** scoring_weights['length_exp']
        
        # Calculate final score
        score = (
            (sequence_coverage ** scoring_weights['coverage_exp']) *
            (avg_tandem_strength ** scoring_weights['tandem_exp']) *
            length_component *
            math.log(total_occurrences + 1) *  # Non-overlapping occurrences
            math.log(max_single_array + 1)
        )
        
        return score
    
    def compile_pattern_statistics(self, pattern, all_data):
        """Compile comprehensive statistics for a pattern"""
        canonical_pattern = self.get_canonical_rotation(pattern)
        
        stats = {
            'pattern': canonical_pattern,
            'pattern_length': len(canonical_pattern),
            'sequences_present': [],
            'sequence_coverage': 0,
            'total_occurrences': 0,  # Non-overlapping
            'total_tandem_arrays': 0,
            'total_tandem_copies': 0,
            'max_tandem_array': 0,
            'avg_tandem_strength': 0,
            'sequence_details': {}
        }
        
        for seq_id, seq_data in all_data.items():
            if canonical_pattern in seq_data:
                data = seq_data[canonical_pattern]
                stats['sequences_present'].append(seq_id)
                stats['total_occurrences'] += data['occurrences']
                stats['total_tandem_arrays'] += data['tandem_array_count']
                stats['total_tandem_copies'] += data['total_tandem_copies']
                stats['max_tandem_array'] = max(stats['max_tandem_array'], data['max_tandem_copies'])
                
                stats['sequence_details'][seq_id] = {
                    'occurrences': data['occurrences'],
                    'tandem_arrays': data['tandem_array_count'],
                    'tandem_copies': data['total_tandem_copies'],
                    'max_array_size': data['max_tandem_copies'],
                    'array_positions': data['tandem_arrays']
                }
        
        stats['sequence_coverage'] = len(stats['sequences_present']) / len(self.sequences)
        stats['avg_tandem_strength'] = stats['total_tandem_copies'] / max(stats['total_tandem_arrays'], 1)
        
        return stats
    
    def find_best_cross_sequence_patterns(self, min_length: int = 100, max_length: int = 200,
                                        min_sequences: int = 2, top_n: int = 20,
                                        scoring_weights=None, step_size: int = 1):
        """Find the best patterns that appear across multiple sequences"""
        print("="*90)
        print("CROSS-SEQUENCE REPEAT PATTERN ANALYSIS")
        print("="*90)
        print(f"Pattern length range: {min_length}-{max_length} bp")
        print(f"Step size: {step_size} bp")
        print(f"Minimum sequences: {min_sequences}")
        
        if scoring_weights:
            print(f"Custom scoring weights: {scoring_weights}")
        
        # Analyze all sequences
        all_data = self.analyze_all_sequences(min_length, max_length, step_size)
        
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
                score = self.score_cross_sequence_pattern(pattern, all_data, scoring_weights)
                stats = self.compile_pattern_statistics(pattern, all_data)
                scored_patterns.append((pattern, stats, score))
        
        print(f"Found {len(scored_patterns)} patterns in {min_sequences}+ sequences")
        
        # Sort by score and return top patterns
        scored_patterns.sort(key=lambda x: x[2], reverse=True)
        return scored_patterns[:top_n]
    
    def visualize_pattern_across_sequences(self, pattern: str, stats: Dict, context: int = 50):
        """Visualize how a pattern appears across sequences"""
        print(f"\n{'='*90}")
        print(f"PATTERN VISUALIZATION: {pattern[:40]}{'...' if len(pattern) > 40 else ''}")
        print(f"{'='*90}")
        print(f"Length: {stats['pattern_length']} bp")
        print(f"Present in: {len(stats['sequences_present'])}/{len(self.sequences)} sequences")
        print(f"Total occurrences (non-overlapping): {stats['total_occurrences']}")
        print(f"Total tandem arrays: {stats['total_tandem_arrays']}")
        print(f"Total tandem copies: {stats['total_tandem_copies']}")
        print(f"Largest array: {stats['max_tandem_array']} copies")
        
        for seq_id in stats['sequences_present']:
            print(f"\n--- {seq_id} ---")
            seq_details = stats['sequence_details'][seq_id]
            print(f"Occurrences: {seq_details['occurrences']}, Arrays: {seq_details['tandem_arrays']}, Max array: {seq_details['max_array_size']} copies")
            
            sequence = self.sequences[seq_id]
            
            for start, end, copies in seq_details['array_positions']:
                array_length = end - start
                print(f"  Array: pos {start:,}-{end:,} ({array_length:,} bp, {copies} copies)")
                
                context_start = max(0, start - context)
                context_end = min(len(sequence), end + context)
                
                before = sequence[context_start:start]
                array_seq = sequence[start:end]
                after = sequence[end:context_end]
                
                print(f"    ...{before[-30:]}[{array_seq[:60]}{'...' if len(array_seq) > 60 else ''}]{after[:30]}...")
    
    def search_for_known_pattern(self, pattern: str, allow_mismatches: int = 0):
        """Search for a specific known pattern across all sequences"""
        print(f"\n{'='*90}")
        print(f"SEARCHING FOR KNOWN PATTERN")
        print(f"{'='*90}")
        print(f"Pattern length: {len(pattern)} bp")
        print(f"Pattern preview: {pattern[:60]}...")
        print(f"Allow mismatches: {allow_mismatches}")
        
        canonical = self.get_canonical_rotation(pattern)
        all_rotations = self.get_all_rotations(pattern)
        
        print(f"Canonical form: {canonical[:60]}...")
        print(f"Total rotations: {len(all_rotations)}")
        
        results = {}
        
        for seq_id, sequence in self.sequences.items():
            print(f"\nSearching in {seq_id} ({len(sequence):,} bp)...")
            
            matches = []
            
            if allow_mismatches == 0:
                # Exact match search for all rotations
                for rotation in all_rotations:
                    pos = 0
                    while True:
                        found = sequence.find(rotation, pos)
                        if found == -1:
                            break
                        matches.append(found)
                        pos = found + 1
            else:
                # Approximate match
                pattern_len = len(pattern)
                for i in range(len(sequence) - pattern_len + 1):
                    candidate = sequence[i:i + pattern_len]
                    
                    for rotation in all_rotations:
                        mismatches = sum(1 for a, b in zip(candidate, rotation) if a != b)
                        if mismatches <= allow_mismatches:
                            matches.append(i)
                            break
            
            matches = sorted(set(matches))
            print(f"  Found {len(matches)} total occurrences")
            
            if matches:
                non_overlapping = self.get_non_overlapping_positions(matches, len(pattern))
                print(f"  Non-overlapping: {len(non_overlapping)}")
                
                tandem_arrays = self.find_tandem_arrays_optimized(matches, len(pattern))
                if tandem_arrays:
                    print(f"  Tandem arrays: {len(tandem_arrays)}")
                    for start, end, copies in tandem_arrays:
                        print(f"    Array at {start:,}-{end:,}: {copies} copies")
                
                results[seq_id] = {
                    'all_positions': matches,
                    'non_overlapping_positions': non_overlapping,
                    'tandem_arrays': tandem_arrays
                }
                
                for i, pos in enumerate(non_overlapping[:3]):
                    context_start = max(0, pos - 30)
                    context_end = min(len(sequence), pos + len(pattern) + 30)
                    before = sequence[context_start:pos]
                    match = sequence[pos:pos + len(pattern)]
                    after = sequence[pos + len(pattern):context_end]
                    print(f"    Match {i+1} at {pos:,}: ...{before[-20:]}[{match[:40]}...]{after[:20]}...")
        
        total_sequences = len([s for s in results if results[s]['all_positions']])
        total_occurrences = sum(len(r['all_positions']) for r in results.values())
        total_non_overlapping = sum(len(r['non_overlapping_positions']) for r in results.values())
        
        print(f"\n{'='*90}")
        print(f"SUMMARY")
        print(f"{'='*90}")
        print(f"Pattern found in: {total_sequences}/{len(self.sequences)} sequences")
        print(f"Total occurrences: {total_occurrences}")
        print(f"Non-overlapping occurrences: {total_non_overlapping}")
        
        return results


def analyze_chromosome_ends(sequences: Dict[str, str], min_pattern: int = 100, max_pattern: int = 200,
                          scoring_weights=None, step_size: int = 1):
    """Main analysis function for chromosome end repeat patterns"""
    print("CHROMOSOME END REPEAT ANALYSIS")
    print("="*90)
    
    finder = ChromosomeEndRepeatFinder(sequences)
    
    results = finder.find_best_cross_sequence_patterns(
        min_length=min_pattern,
        max_length=max_pattern,
        min_sequences=2,
        top_n=15,
        scoring_weights=scoring_weights,
        step_size=step_size
    )
    
    if not results:
        print("No cross-sequence patterns found!")
        return finder, None, results
    
    print(f"\n{'='*90}")
    print("TOP CROSS-SEQUENCE PATTERNS")
    print("="*90)
    print(f"{'Rank':>4} {'Length':>6} {'Seqs':>4} {'Occur':>6} {'Arrays':>6} {'Copies':>7} {'Max':>5} {'Score':>10}")
    print("-"*90)
    
    for i, (pattern, stats, score) in enumerate(results):
        print(f"{i+1:4d} {stats['pattern_length']:6d} {len(stats['sequences_present']):4d} "
              f"{stats['total_occurrences']:6d} {stats['total_tandem_arrays']:6d} "
              f"{stats['total_tandem_copies']:7d} {stats['max_tandem_array']:5d} {score:10.2f}")
    
    best_pattern, best_stats, best_score = results[0]
    finder.visualize_pattern_across_sequences(best_pattern, best_stats)
    
    return finder, best_pattern, results


# Test function
def create_test_chromosome_sequences():
    """Create test sequences with large patterns that form tandem repeats"""
    large_pattern = (
        "TTAGGGTTACGTACGTACGTACGTACGTACGTACGTACGT"
        "AAAAAGGGGCCCCTTTTTATCGATCGATCGATCGATCGAT"
        "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
        "NNNATCGATCG"
    )
    
    sequences = {}
    
    seq1 = "N" * 200
    seq1 += large_pattern * 5
    seq1 += "ATCGATCGATCG" * 20
    seq1 += large_pattern * 3
    seq1 += "N" * 150
    sequences["chr1p_telomere"] = seq1
    
    seq2 = "N" * 180
    seq2 += large_pattern * 7
    seq2 += "GCTAGCTA" * 25
    seq2 += large_pattern * 2
    seq2 += "N" * 100
    sequences["chr2q_telomere"] = seq2
    
    rotated_pattern = large_pattern[50:] + large_pattern[:50]
    seq3 = "N" * 220
    seq3 += rotated_pattern * 4
    seq3 += "AAAATTTT" * 30
    seq3 += rotated_pattern * 6
    seq3 += "N" * 80
    sequences["chr3p_telomere"] = seq3
    
    return sequences, large_pattern


def execute(telomeres: List[Optional[TelomereSequence]]) -> Optional[str]:
    test_dict: Dict[str, str] = {}
    for telomer in telomeres:
        if telomer and telomer.sequence:
            test_dict[telomer.chromosome_end_id] = telomer.sequence

    finder, best_pattern, results = analyze_chromosome_ends(test_dict, 160, 170)


    print(f"best_pattern: {best_pattern}")
    for r in results:
        print(r)
    

    if best_pattern:
      return best_pattern
    return None


if __name__ == "__main__":
    print("Creating test sequences...")
    test_sequences, expected_pattern = create_test_chromosome_sequences()
    
    print(f"\nExpected pattern length: {len(expected_pattern)} bp")
    
    finder, best_pattern, results = analyze_chromosome_ends(test_sequences, 100, 180, step_size=1)
    
    if results:
        print(f"\n{'='*90}")
        print("VALIDATION")
        print(f"{'='*90}")
        print(f"Expected pattern: {len(expected_pattern)} bp")
        print(f"Best pattern found: {len(results[0][0])} bp")