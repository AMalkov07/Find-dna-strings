from utils.data_structures import TelomereSequence
from typing import List

class AlignmentStrategy:
    def __init__(self, telomeres: List[TelomereSequence], pattern: str):
        self.telomers = telomeres
        self. pattern = pattern
    
    def execute(self):
        print("executing alignment strategy")