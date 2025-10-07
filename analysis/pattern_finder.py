from typing import List, Optional, Dict, Tuple
from collections import defaultdict
import math
import sys

from utils.data_structures import Config, TelomereSequence


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
    
    def find_best_repeats_efficient(self, min_length: int = 300, max_length: int = 400, 
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
        for c in candidates:
            print(c)
        sys.exit()
        return candidates



class PatternFinder:
    def __init__(self, telomeres: List[Optional[TelomereSequence]], config: Config):
        self.telomers: List[Optional[TelomereSequence]] = telomeres
        self.config = config


    def execute(self) -> str:
        for telomer in self.telomers:
            if telomer and telomer.sequence:
                finder = EfficientDNARepeatFinder(telomer.sequence)
                result = finder.find_best_repeats_efficient()
                
        return "done"