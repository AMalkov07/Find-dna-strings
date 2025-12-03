
from utils.data_structures import TelomereSequence, Config, TemplateSwitchEvent, TemplateSwitchData
from typing import List, Optional, Tuple


class MatchResult:
    """Represents the result of finding a position in the circular DNA."""
    
    def __init__(self, match_type, positions, length):
        self.type = match_type  # 'unique', 'ambiguous', or 'mutation'
        self.positions = positions  # List of circle positions (empty for mutations)
        self.length = length
    
    def __repr__(self):
        return f"MatchResult(type={self.type}, positions={self.positions}, length={self.length})"


class TemplateSwitchingStrategy:
    """
    Analyzes a long DNA string against a circular reference to identify
    uninterrupted sections without jumps in the circular sequence.
    """
    
    def __init__(self, telomeres: List[Optional[TelomereSequence]], circular_dna: str, config: Config):
        """
        Initialize the analyzer.
        
        Args:
            circular_dna: The circular reference DNA string (~150 bp)
            long_dna: The longer DNA string composed of sections from circular_dna
        """
        self.telomers: List[Optional[TelomereSequence]] = telomeres
        self.circular = circular_dna
        self.circle_len = len(circular_dna)
        self.config = config
        self.current_analysis = TemplateSwitchData([], [])
        self.mutation_lookahead = config.mutation_lookahead
    
    def find_uninterrupted_segments(self):
        """
        Identify all uninterrupted segments in the long DNA string.
        
        Returns:
            List of Segment objects
        """

        added_first = False

        self.segments = []
        i = 0
        
        while i < len(self.long):
            # Find the position and type for current position in long string
            result = self._find_unique_position(i)
            
            if result is None:
                i += 1
                continue
            
            if result.type == 'unique':
                # Unique position found
                event = self._extend_match(i, result.positions[0], result.length)
            elif result.type == 'ambiguous':
                # Non-unique but matchable position (ambiguous)
                # We'll extend from one of the possible positions to get the length
                event = self._extend_match(i, result.positions[0], result.length, 
                                            ambiguous=True, possible_positions=result.positions)
            else: 
                event = TemplateSwitchEvent(self.long[i], i, i, None, None, True, [], [], [])
            
            if added_first or len(event.telomer_chunk) >= min(50, self.circle_len//2):
                if not added_first:
                    added_first = True
                self.current_analysis.template_switch_event_indexes.append(i)
                self.current_analysis.template_switch_events.append(event)
            i = event.telomer_end+1 #make sure that the +! is correct***********
        
        return self.segments
    
    def _find_unique_position(self, long_idx):
        """
        Find the minimum length substring starting at long_idx that maps
        to position(s) in the circular DNA.
        
        Returns:
            MatchResult object or None
        """
        if long_idx >= len(self.long):
            return None
        
        ambiguous_result = None
        
        # Start with length 1 and increase until we find matches
        for test_len in range(1, min(len(self.long) - long_idx + 1, self.circle_len + 1)):
            substring = self.long[long_idx:long_idx + test_len]
            matches = self._find_all_matches_in_circle(substring)
            
            if len(matches) == 1:
                # Found unique match
                return MatchResult('unique', matches, test_len)
            elif len(matches) > 1:
                # Multiple matches - continue to see if we can make it unique
                # But remember this as our fallback
                ambiguous_result = MatchResult('ambiguous', matches, test_len)
            elif len(matches) == 0:
                # No match found
                if test_len == 1:
                    # Single base pair not in circle = mutation
                    return MatchResult('mutation', [], 1)
                else:
                    # Was ambiguous before, now no matches - return previous ambiguous
                    if ambiguous_result is not None:
                        return ambiguous_result
                    else:
                        # Length > 1 but never had matches - shouldn't happen normally
                        # but treat first bp as mutation
                        return MatchResult('mutation', [], 1)
        
        # If we get here, even the maximum length isn't unique
        # Return the ambiguous result if we found one
        if ambiguous_result is not None:
            return ambiguous_result
        
        # Shouldn't reach here, but handle gracefully
        return MatchResult('mutation', [], 1)
    
    def _find_all_matches_in_circle(self, substring) -> List[int]:
        """
        Find all positions in the circular DNA where substring matches.
        Handles wrapping around the circle.
        
        Returns:
            List of starting positions in circular DNA
        """
        matches = []
        sub_len = len(substring)
        
        # Create a doubled circular string to handle wrapping
        doubled_circle = self.circular + self.circular
        
        for i in range(self.circle_len):
            if doubled_circle[i:i + sub_len] == substring:
                matches.append(i)
        
        # If there are exactly 2 matches, they could be the same position
        # if the distance between the positions is equal to the circle size, then they should be the same position
        if len(matches) == 2:
            if abs(matches[1] - matches[0]) == self.circle_len:
                return [matches[0]]
        
        return matches

        
    def _detect_mutation(self, long_pos, circle_pos):
        """
        Detect if there's a recoverable mutation (mismatch, insertion, or deletion).
        
        Args:
            long_pos: Current position in long DNA
            circle_pos: Current position in circular DNA
        
        Returns:
            Dictionary with mutation info or None if no recoverable mutation
        """
        lookahead = self.mutation_lookahead
        
        # Make sure we have enough room to look ahead
        if long_pos + lookahead > len(self.long):
            return None
        
        # Check for single base mismatch
        # long[i+1] != circle[j+1], but long[i+2:i+11] == circle[j+2:j+11]
        mismatch_match = True
        for k in range(1, lookahead):
            if long_pos + k + 1 >= len(self.long):
                mismatch_match = False
                break
            if self.long[long_pos + k + 1] != self.circular[(circle_pos + k + 1) % self.circle_len]:
                mismatch_match = False
                break
        
        if mismatch_match:
            return {
                'type': 'mismatch',
                'long_pos': long_pos,
                'circle_pos': circle_pos,
                'long_base': self.long[long_pos],
                'circle_base': self.circular[circle_pos % self.circle_len]
            }
        
        # Check for insertion
        # long[i+1] != circle[j+1], but long[i+2:i+11] == circle[j+1:j+10]
        insertion_match = True
        for k in range(lookahead - 1):
            if long_pos + k + 2 >= len(self.long):
                insertion_match = False
                break
            if self.long[long_pos + k + 2] != self.circular[(circle_pos + k + 1) % self.circle_len]:
                insertion_match = False
                break
        
        if insertion_match:
            return {
                'type': 'insertion',
                'long_pos': long_pos,
                'circle_pos': circle_pos,
                'inserted_base': self.long[long_pos]
            }
        
        # Check for deletion
        # long[i+1] != circle[j+1], but long[i+1:i+10] == circle[j+2:j+11]
        deletion_match = True
        for k in range(lookahead - 1):
            if long_pos + k + 1 >= len(self.long):
                deletion_match = False
                break
            if self.long[long_pos + k + 1] != self.circular[(circle_pos + k + 2) % self.circle_len]:
                deletion_match = False
                break
        
        if deletion_match:
            return {
                'type': 'deletion',
                'long_pos': long_pos,
                'circle_pos': circle_pos,
                'deleted_base': self.circular[circle_pos % self.circle_len]
            }
        
        return None


    #def _extend_match(self, long_start, circle_start, min_length, 
                     #ambiguous=False, possible_positions=None):
        #"""
        #Extend a match as far as possible without jumps.
        #We already know the first min_length characters match uniquely.
        
        #Returns a Segment object.
        #"""
        #length = min_length
        #circle_pos = circle_start + min_length
        #long_pos = long_start + min_length
        #mutations = []
        
        ## Extend while characters match and follow circular sequence
        #while long_pos < len(self.long):
            #if self.long[long_pos] != self.circular[circle_pos % self.circle_len]:
                #break
            
            #length += 1
            #long_pos += 1
            #circle_pos += 1
        
        ## Get the sequence
        #sequence = self.long[long_start:long_start + length]
        
        #if ambiguous:
            #return TemplateSwitchEvent(
                #sequence,
                #long_start,
                #long_start + length,
                #'ambiguous',
                #'ambiguous',
                #False
            #)
        #else:
            ## Determine actual circle end position and if it wraps
            #actual_circle_end = circle_start + length
            #wraps = actual_circle_end >= self.circle_len
            #circle_end = actual_circle_end % self.circle_len
            
            ## If circle_end is 0 and we wrapped, it means we ended exactly at the circle length
            #if circle_end == 0 and wraps:
                #circle_end = self.circle_len
            
            #return TemplateSwitchEvent(
                #sequence,
                #long_start,
                #long_start + length,
                #circle_start,
                #circle_end,
                #False
            #)
    
    def _extend_match(self, long_start, circle_start, min_length, 
                     ambiguous=False, possible_positions=None):
        """
        Extend a match as far as possible, detecting mutations along the way.
        We already know the first min_length characters match.
        
        Returns a Segment object.
        """
        length = min_length
        circle_pos = circle_start + min_length
        long_pos = long_start + min_length
        insertion_events: List[Tuple[int, str]] = []
        deletion_events: List[Tuple[int, str]] = []
        mismatch_events: List[Tuple[int, str, str]] = []
        
        # Extend while characters match (or we can detect mutations)
        while long_pos < len(self.long):
            long_base = self.long[long_pos]
            circle_base = self.circular[circle_pos % self.circle_len]
            
            if long_base == circle_base:
                # Perfect match, continue
                length += 1
                long_pos += 1
                circle_pos += 1
            else:
                # Mismatch detected - check for mutation types
                mutation = self._detect_mutation(long_pos, circle_pos % self.circle_len)
                
                if mutation is None:
                    # No recoverable mutation, end of segment
                    break
                
                if mutation['type'] == 'mismatch':
                    # Single base mismatch: both advance by 1
                    length += 1
                    long_pos += 1
                    circle_pos += 1
                    mismatch_events.append((mutation['long_pos'], mutation['long_base'], mutation['circle_base']))
                elif mutation['type'] == 'insertion':
                    # Insertion in long: long advances, circle stays
                    length += 1
                    long_pos += 1
                    insertion_events.append((mutation['long_pos'], mutation['inserted_base']))
                    # circle_pos stays the same
                elif mutation['type'] == 'deletion':
                    # Deletion in long: circle advances, long stays
                    # Don't increment length (no new base in long)
                    circle_pos += 1
                    deletion_events.append((mutation['circle_pos'], mutation['deleted_base']))
                    # long_pos stays the same
        
        # Get the sequence
        sequence = self.long[long_start:long_start + length]
        
        if ambiguous:
            return TemplateSwitchEvent(
                sequence,
                long_start,
                long_start + length,
                'ambiguous',
                'ambiguous',
                False,
                insertion_events,
                deletion_events,
                mismatch_events
            )
        else:
            # Determine actual circle end position and if it wraps
            actual_circle_end = circle_pos
            wraps = actual_circle_end >= self.circle_len
            circle_end = actual_circle_end % self.circle_len
            
            # If circle_end is 0 and we wrapped, it means we ended exactly at the circle length
            if circle_end == 0 and wraps:
                circle_end = self.circle_len
            
            return TemplateSwitchEvent(
                sequence,
                long_start,
                long_start + length,
                circle_start,
                circle_end,
                False,
                insertion_events,
                deletion_events,
                mismatch_events
            )
    
    def print_segments(self):
        """Print a formatted summary of all segments."""
        if not self.segments:
            print("No segments found. Run find_uninterrupted_segments() first.")
            return
        
        print(f"Found {len(self.segments)} segments:\n")
        
        for idx, seg in enumerate(self.segments, 1):
            print(f"Segment {idx} ({seg.type.upper()}):")
            print(f"  Long DNA: [{seg.long_start}:{seg.long_end}]")
            
            if seg.type == 'unique':
                print(f"  Circle DNA: [{seg.circle_start}:{seg.circle_end}]")
                print(f"  Wraps around circle: {seg.wraps}")
            elif seg.type == 'ambiguous':
                print(f"  Circle DNA: [ambiguous]")
                print(f"  Possible positions: {seg.possible_positions}")
            else:  # mutation
                print(f"  Circle DNA: [mutation - not found in circle]")
            
            print(f"  Length: {seg.length} bp")
            print(f"  Sequence: {seg.sequence[:50]}{'...' if len(seg.sequence) > 50 else ''}")
            print()
    
    def get_coverage_stats(self):
        """Calculate statistics about circle coverage."""
        if not self.segments:
            return None
        
        total_length = sum(seg.length for seg in self.segments)
        long_coverage = (total_length / len(self.long)) * 100 if len(self.long) > 0 else 0
        
        num_unique = sum(1 for seg in self.segments if seg.type == 'unique')
        num_ambiguous = sum(1 for seg in self.segments if seg.type == 'ambiguous')
        num_mutations = sum(1 for seg in self.segments if seg.type == 'mutation')
        
        return {
            'num_segments': len(self.segments),
            'num_unique': num_unique,
            'num_ambiguous': num_ambiguous,
            'num_mutations': num_mutations,
            'total_mapped_length': total_length,
            'long_dna_coverage': long_coverage,
            'avg_segment_length': total_length / len(self.segments)
        }


    def execute(self) -> None:
        #analysis_output: List[Optional[TemplateSwitchData]] = [None] * self.config.max_ends
        for telomer in self.telomers:
            if telomer and telomer.sequence:
                self.long = telomer.sequence
                self.find_uninterrupted_segments()
                telomer.analysis = self.current_analysis
                self.current_analysis = TemplateSwitchData([], [])
        #return analysis_output