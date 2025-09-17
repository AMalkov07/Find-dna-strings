from dataclasses import dataclass
from typing import List, Dict, Tuple, Any, Optional

@dataclass
class TelomereSequence:
    survivor_name: str
    sequence: Optional[str]
    chromosome_end_id: str
    
@dataclass
class AnalysisResult:
    telomere_id: str
    position: int
    substring: str
    result_type: str  # "alignment" or "template_switch"
    details: Dict[str, Any]  # alignment scores, switch points, etc.

@dataclass
class Config:
    fasta_file_path: str
    output_file: str
    analysis_strategy: str
    max_ends: int
    pattern: str
    maximum_alignment_mutations: int
    skip_seeding: Optional[bool]

@dataclass
class TemplateSwitchEvent:
    telomer_chunk: str
    telomer_start: int
    telomer_end: int
    pattern_start: Any
    pattern_end: Any
    is_mutation: bool


#fix: find a beetter way to connect the 2 arrays
@dataclass
class TemplateSwitchData:
    template_switch_event_indexes: List[int]
    template_switch_events: List[TemplateSwitchEvent]

@dataclass
class ImperfectAlignmentEvent:
    full_telomer_start_index: Optional[List[int]]
    mutagenic_zone_start_index: int
    mutagenic_zone_end_index: int
    insertion_events: List[Tuple[int, str]]
    deletion_events: List[Tuple[int, str]]
    mismatch_events: List[Tuple[int, str, str]]
    score: float

@dataclass
class AlignmentData:
    perfect_alignments: List[int]
    imperfect_alignments: List[ImperfectAlignmentEvent]

@dataclass
class AlignmentSettings:
    match: int
    mismatch: int
    gap_open: int
    gap_extend: int

@dataclass
class SeedExtendSettings:
    alignment_settings: AlignmentSettings
    k: int
    flank: int
    min_identity: float
    offset_tolerance: int
    min_seed: int
    last_zone: bool

@dataclass
class SeedExtendCluster:
    offset: int
    hits: List[Tuple[int, int]]
    qmin: int
    qmax: int
    rmin: int
    rmax: int
    n_seeds: int