from dataclasses import dataclass
from typing import List, Dict, Any

@dataclass
class TelomereSequence:
    survivor_name: str
    sequence: str
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

@dataclass
class TemplateSwitchEvent:
    telomer_chunk: str
    telomer_start: int
    telomer_end: int
    pattern_start: int
    pattern_end: int
    is_mutation: bool


#fix: find a beetter way to connect the 2 arrays
@dataclass
class TemplateSwitchData:
    template_switch_event_indexes: List[int]
    template_switch_events: List[TemplateSwitchEvent]