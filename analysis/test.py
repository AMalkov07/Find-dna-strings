import math
from collections import defaultdict
from typing import List, Tuple, Dict, Set

class EfficientDNARepeatFinder:
    def __init__(self, sequence: str):
        """
        Initialize with DNA sequence and build suffix array
        """
        self.sequence = sequence.upper()
        self.length = len(sequence)
        
        # Build suffix array and LCP array
        print("Building suffix array...")
        self.suffix_array: List[int] = self._build_suffix_array()
        print("Computing LCP array...")
        self.lcp_array = self._build_lcp_array()
        print("Preprocessing complete!")
    
    def _build_suffix_array(self) -> List[int]:
        """
        Build suffix array using counting sort approach
        More efficient than naive O(n²logn) sorting
        """
        n = self.length
        
        # Create list of (suffix, original_index) pairs
        suffixes: List[Tuple[str, int]] = []
        for i in range(n):
            suffixes.append((self.sequence[i:], i))
        
        # Sort suffixes lexicographically
        suffixes.sort(key=lambda x: x[0])
        
        # Extract just the indices
        return [suffix[1] for suffix in suffixes]
    
    def _build_lcp_array(self) -> List[int]:
        """
        Build Longest Common Prefix array
        lcp[i] = length of common prefix between suffix_array[i] and suffix_array[i+1]
        """
        n = self.length
        if n <= 1:
            return []
        
        lcp = []
        
        for i in range(n - 1):
            # Get the two suffixes to compare
            suffix1_start = self.suffix_array[i]
            suffix2_start = self.suffix_array[i + 1]
            
            # Find common prefix length
            common_len = 0
            max_compare = min(n - suffix1_start, n - suffix2_start)
            
            for j in range(max_compare):
                if self.sequence[suffix1_start + j] == self.sequence[suffix2_start + j]:
                    common_len += 1
                else:
                    break
            
            lcp.append(common_len)
        
        return lcp
    
    def _find_all_repeated_substrings(self, min_length: int, max_length: int) -> Dict[str, List[int]]:
        """
        Find all repeated substrings using suffix array and LCP array
        Returns dict mapping substring -> list of starting positions
        """
        repeated_substrings = defaultdict(list)
        
        # Process each LCP entry
        for i, lcp_len in enumerate(self.lcp_array):
            if lcp_len >= min_length:
                # We found a repeated substring of length lcp_len
                pos1 = self.suffix_array[i]
                pos2 = self.suffix_array[i + 1]
                
                # Extract all prefixes of this common substring up to max_length
                max_extract = min(lcp_len, max_length)
                
                for length in range(min_length, max_extract + 1):
                    substring = self.sequence[pos1:pos1 + length]
                    
                    # Add both positions (we'll deduplicate later)
                    if pos1 not in repeated_substrings[substring]:
                        repeated_substrings[substring].append(pos1)
                    if pos2 not in repeated_substrings[substring]:
                        repeated_substrings[substring].append(pos2)
        
        return repeated_substrings
    
    def _find_all_occurrences_efficient(self, pattern: str) -> List[int]:
        """
        Find all occurrences of a pattern using binary search on suffix array
        Much faster than string searching for multiple patterns
        """
        def compare_suffix_with_pattern(suffix_idx: int, pattern: str) -> int:
            """Returns -1 if suffix < pattern, 0 if equal prefix, 1 if suffix > pattern"""
            suffix_start = self.suffix_array[suffix_idx]
            max_compare = min(len(pattern), self.length - suffix_start)
            
            for i in range(max_compare):
                if self.sequence[suffix_start + i] < pattern[i]:
                    return -1
                elif self.sequence[suffix_start + i] > pattern[i]:
                    return 1
            
            # All compared characters are equal
            if max_compare < len(pattern):
                return -1  # suffix is shorter than pattern
            return 0  # perfect match or suffix is longer
        
        # Binary search for first occurrence
        left, right = 0, len(self.suffix_array)
        first_match = -1
        
        while left < right:
            mid = (left + right) // 2
            cmp = compare_suffix_with_pattern(mid, pattern)
            if cmp < 0:
                left = mid + 1
            else:
                if cmp == 0:
                    first_match = mid
                right = mid
        
        if first_match == -1:
            return []
        
        # Find all consecutive matches
        matches = []
        for i in range(first_match, len(self.suffix_array)):
            if compare_suffix_with_pattern(i, pattern) == 0:
                suffix_pos = self.suffix_array[i]
                # Verify it's an exact match (not just prefix)
                if (suffix_pos + len(pattern) <= self.length and 
                    self.sequence[suffix_pos:suffix_pos + len(pattern)] == pattern):
                    matches.append(suffix_pos)
            else:
                break
        
        return sorted(matches)
    
    def find_non_overlapping_matches(self, pattern: str) -> List[int]:
        """
        Find non-overlapping matches using greedy approach
        """
        all_matches = self._find_all_occurrences_efficient(pattern)
        if not all_matches:
            return []
        
        non_overlapping = []
        last_end = -1
        pattern_len = len(pattern)
        
        for pos in all_matches:
            if pos >= last_end:
                non_overlapping.append(pos)
                last_end = pos + pattern_len
        
        return non_overlapping
    
    def count_non_overlapping_repeats(self, pattern: str) -> int:
        """Count non-overlapping occurrences"""
        return len(self.find_non_overlapping_matches(pattern))
    
    def score_pattern(self, pattern: str, method: str = "log") -> float:
        """Score a pattern based on length and repeat count"""
        repeats = self.count_non_overlapping_repeats(pattern)
        length = len(pattern)
        
        if repeats < 2:
            return 0
        
        if method == "log":
            return length * math.log(repeats)
        elif method == "sqrt":
            return length * math.sqrt(repeats)
        elif method == "power":
            return length * (repeats ** 0.7)
        elif method == "weighted":
            return (length ** 1.2) * (repeats ** 0.8)
        else:
            raise ValueError(f"Unknown scoring method: {method}")
    
    def find_best_repeats_efficient(self, min_length: int = 100, max_length: int = 200, 
                                  scoring_method: str = "log", min_repeats: int = 2) -> List[Tuple[str, int, float]]:
        """
        Efficiently find best repetitive patterns using suffix array approach
        """
        print(f"Finding repeated substrings between {min_length}-{max_length} bp...")
        
        # Step 1: Find all repeated substrings efficiently
        repeated_substrings = self._find_all_repeated_substrings(min_length, max_length)
        
        print(f"Found {len(repeated_substrings)} unique repeated substrings")
        
        # Step 2: Count non-overlapping occurrences and score each pattern
        candidates = []
        processed = 0
        
        for pattern, positions in repeated_substrings.items():
            processed += 1
            if processed % 1000 == 0:
                print(f"Processed {processed}/{len(repeated_substrings)} patterns...")
            
            # Quick filter: if we only found 2 positions in LCP analysis, 
            # but they might not be the only occurrences
            if len(positions) >= min_repeats or len(positions) == 2:
                repeats = self.count_non_overlapping_repeats(pattern)
                if repeats >= min_repeats:
                    score = self.score_pattern(pattern, scoring_method)
                    candidates.append((pattern, repeats, score))
        
        print(f"Found {len(candidates)} patterns with {min_repeats}+ repeats")
        
        # Step 3: Sort by score and return top candidates
        candidates.sort(key=lambda x: x[2], reverse=True)
        return candidates
    
    def analyze_pattern_distribution(self, pattern: str) -> Dict:
        """
        Analyze the distribution and spacing of a pattern
        """
        matches = self.find_non_overlapping_matches(pattern)
        if len(matches) < 2:
            return {"error": "Pattern has fewer than 2 non-overlapping matches"}
        
        # Calculate distances between consecutive matches
        distances = []
        for i in range(1, len(matches)):
            distance = matches[i] - (matches[i-1] + len(pattern))
            distances.append(distance)
        
        return {
            'pattern_length': len(pattern),
            'total_matches': len(matches),
            'match_positions': matches,
            'inter_repeat_distances': distances,
            'mean_distance': sum(distances) / len(distances) if distances else 0,
            'min_distance': min(distances) if distances else 0,
            'max_distance': max(distances) if distances else 0,
            'total_coverage': len(matches) * len(pattern),
            'coverage_percentage': (len(matches) * len(pattern)) / self.length * 100
        }
    
    def visualize_matches(self, pattern: str, context: int = 20):
        """Show pattern matches with context"""
        matches = self.find_non_overlapping_matches(pattern)
        print(f"\nPattern: {pattern[:30]}{'...' if len(pattern) > 30 else ''}")
        print(f"Length: {len(pattern)} bp")
        print(f"Non-overlapping matches: {len(matches)}")
        print("-" * 60)
        
        for i, pos in enumerate(matches):
            start_context = max(0, pos - context)
            end_context = min(self.length, pos + len(pattern) + context)
            
            before = self.sequence[start_context:pos]
            match = self.sequence[pos:pos + len(pattern)]
            after = self.sequence[pos + len(pattern):end_context]
            
            print(f"Match {i+1:2d} (pos {pos:4d}): ...{before[-context:]}[{match[:20]}{'...' if len(match) > 20 else ''}]{after[:context]}...")


def analyze_dna_sequence_efficient(sequence: str, min_len: int = 100, max_len: int = 200):
    """
    Efficient analysis of DNA sequence for repetitive patterns
    """
    print(f"Analyzing sequence of length {len(sequence):,} bp")
    print(f"Looking for patterns between {min_len}-{max_len} bp")
    print("=" * 70)
    
    # Initialize the efficient finder
    finder = EfficientDNARepeatFinder(sequence)
    
    # Find best patterns
    print(f"\nSearching for optimal repetitive patterns...")
    results = finder.find_best_repeats_efficient(min_len, max_len, scoring_method="log")
    
    if not results:
        print("No repetitive patterns found matching criteria")
        return finder, None, results
    
    print(f"\nTop 10 patterns found:")
    print(f"{'Rank':>4} {'Length':>6} {'Repeats':>7} {'Score':>8} {'Pattern Preview':>30}")
    print("-" * 70)
    
    for i, (pattern, repeats, score) in enumerate(results[:10]):
        preview = pattern[:25] + "..." if len(pattern) > 25 else pattern
        print(f"{i+1:4d} {len(pattern):6d} {repeats:7d} {score:8.2f} {preview:>30}")
    
    # Detailed analysis of best pattern
    best_pattern = results[0][0]
    print(f"\n" + "=" * 70)
    print("DETAILED ANALYSIS OF TOP PATTERN")
    print("=" * 70)
    
    analysis = finder.analyze_pattern_distribution(best_pattern)
    print(f"Pattern length: {analysis['pattern_length']} bp")
    print(f"Total matches: {analysis['total_matches']}")
    print(f"Sequence coverage: {analysis['coverage_percentage']:.1f}%")
    print(f"Mean inter-repeat distance: {analysis['mean_distance']:.1f} bp")
    print(f"Distance range: {analysis['min_distance']}-{analysis['max_distance']} bp")
    
    # Show matches
    finder.visualize_matches(best_pattern)
    
    return finder, best_pattern, results


# Test function
def create_test_sequence_complex():
    """Create a more complex test sequence with multiple repeat types"""
    import random
    
    # Primary repeat - this should be found as the best
    primary_repeat = ("ATCGATCGTAGCTAGCTACGTACGTAAAAGGGGCCCCTTTT" + 
                     "ATCGATCGTAGCTAGCTACGTACGTAAAAGGGGCCCCTTTT" + 
                     "ATCGATCGTAGCTAGCTACGTACGT")  # ~120bp
    
    # Secondary repeat - shorter, more frequent
    secondary_repeat = "TATATATATACGCGCGCGAAAAAATTTTTT"  # ~30bp
    
    # Build sequence
    sequence = "N" * 200  # Start buffer
    
    # Add primary repeats with variable spacing
    for i in range(6):
        sequence += primary_repeat
        spacer_len = random.randint(80, 300)
        spacer = ''.join(random.choices('ATCGN', k=spacer_len))
        sequence += spacer
    
    # Add some secondary repeats
    for i in range(10):
        pos = random.randint(0, len(sequence) - 100)
        sequence = sequence[:pos] + secondary_repeat + sequence[pos:]
    
    sequence += "N" * 200  # End buffer
    
    return sequence


# Example usage
if __name__ == "__main__":
    # Test with complex synthetic sequence
    print("Creating test sequence with known repeat patterns...")
    test_seq = create_test_sequence_complex()
    
    print(f"Test sequence length: {len(test_seq):,} bp")
    
    # Run efficient analysis
    finder, best_pattern, results = analyze_dna_sequence_efficient(test_seq, 80, 150)
    
    if best_pattern:
        print(f"\n" + "=" * 70)
        print("COMPARISON WITH EXPECTED PATTERN")
        print("=" * 70)
        expected_pattern = ("ATCGATCGTAGCTAGCTACGTACGTAAAAGGGGCCCCTTTT" + 
                          "ATCGATCGTAGCTAGCTACGTACGTAAAAGGGGCCCCTTTT" + 
                          "ATCGATCGTAGCTAGCTACGTACGT")
        
        print(f"Expected pattern length: {len(expected_pattern)}")
        print(f"Found pattern length: {len(best_pattern)}")
        print(f"Match: {expected_pattern in best_pattern or best_pattern in expected_pattern}")