import sys
from typing import List, Tuple, Optional
from collections import deque

from utils.data_structures import TelomereSequence, AlignmentData, Config, ImperfectAlignmentEvent, AlignmentSettings, SeedExtendSettings
from analysis.alignment_seed_and_extend_pipeline import SeedAndExtend

class AlignmentStrategy:
    def __init__(self, telomeres: List[TelomereSequence], pattern: str, config: Config):
        self.telomers = telomeres
        self.pattern = pattern
        self.config = config
        self.alignment_settings = AlignmentSettings(match=5, mismatch=-4, gap_open=10, gap_extend=2)

    def _identify_perfect_alignments(self, telomer_str: str) -> List[int]:
        pattern = self.pattern
        n_pattern = len(pattern)
        positions: List[int] = []
        i = 0 #might need to set i to an offset instead
        while i < len(telomer_str) - n_pattern:
            if telomer_str[i:i+n_pattern] == pattern:
                positions.append(i)
                i += n_pattern
            else:
                i += 1
        return positions

    def _get_mutagenic_zones(self, telomer_str: str, perfect_alignments: List[int]) -> Tuple[List[str], List[int]]:
        n_pattern = len(self.pattern)
        str_locations: List[int] = []
        mutagenic_zone: List[str]= []

        if len(perfect_alignments) < 1:
            str_locations.append(0)
            mutagenic_zone.append(telomer_str)
        else:
            last_val = 0
            for i in range(0, len(perfect_alignments)):
                alignment = perfect_alignments[i]
                tmp = telomer_str[last_val:alignment]
                str_locations.append(last_val)
                mutagenic_zone.append(tmp)
                last_val = alignment+n_pattern
            tmp = telomer_str[last_val:]
            str_locations.append(last_val)
            mutagenic_zone.append(tmp)

        return mutagenic_zone, str_locations


    def _identiy_imperfect_alignments(self, telomer_str: str, perfect_alignments: List[int]) -> List[ImperfectAlignmentEvent]:
        n_pattern = len(self.pattern)
        max_mistakes = self.config.maximum_alignment_mutations

        mutagenic_zone, str_locations = self._get_mutagenic_zones(telomer_str, perfect_alignments)
        queue = deque(mutagenic_zone)

        last_elem = False
        counter = -1
        min_mutagenic_zone_len = n_pattern - n_pattern/ max_mistakes

        while queue:
            counter += 1
            current_mutagenic_zone = queue.popleft()
            if len(queue) == 0:
                last_elem = True
            if len(current_mutagenic_zone) < min_mutagenic_zone_len:
                continue

            seed_hit_settings = SeedExtendSettings(alignment_settings=self.alignment_settings, k=15, flank=80, min_identity=.9, offset_tolerance=6, min_seed=2)
            seed_and_extend_pipeline = SeedAndExtend(seed_hit_settings)
            seed_and_extend_pipeline.execute()

    def execute(self) -> List[Optional[AlignmentData]]:
        analysis_output: List[Optional[AlignmentData]] = [None] * self.config.max_ends
        for i, telomer in enumerate(self.telomers):
            if telomer and telomer.sequence:
                perfect_alignments: List[int] = self._identify_perfect_alignments(telomer.sequence)
                analysis_output[i] = AlignmentData(perfect_alignments=perfect_alignments, imperfect_alignments=[])
                imperfect_alignment_indexes, alignment_insertions_and_deletions = self._identiy_imperfect_alignments(telomer.sequence, perfect_alignments)
        return analysis_output