from utils.data_structures import TelomereSequence, Config, TemplateSwitchEvent, TemplateSwitchData
from typing import List, Optional

class TemplateSwitchingStrategy:
    def __init__(self, telomeres: List[Optional[TelomereSequence]], pattern: str, config: Config):
        self.telomers: List[Optional[TelomereSequence]] = telomeres
        self.pattern: str = pattern
        self.doubled_pattern = pattern + pattern
        self.config = config

    def _find_all_starts(self, window: str, offset: int) -> List[int]:
        tmp_doubled_pattern = self.doubled_pattern
        n_pattern = len(self.pattern)
        output_starting_pos: List[int] = []
        while True:
            start = tmp_doubled_pattern.find(window) + 1
            if start == 0:
                break
            start += offset
            output_starting_pos.append(start % n_pattern)
            tmp_doubled_pattern = tmp_doubled_pattern[start + len(window):]
        return output_starting_pos
            
        
    #Fix: read through this again and break up into other functions
    def _identify_template_switches(self, telomer_str: str) -> Optional[TemplateSwitchData]:
        pattern = self.pattern
        doubled_pattern = self.doubled_pattern
        n_telomer = len(telomer_str)
        n_pattern = len(self.pattern)
        offset = 0
        prev_window = ""
        current_analysis = TemplateSwitchData([], [])
        first_unit_added = False
        for i, c in enumerate(telomer_str):
            if c == 'A' and len(prev_window) == 0:
                #self.my_dict[key].template_switching_indexes.append(i + offset)
                if first_unit_added:
                    current_analysis.template_switch_event_indexes.append(i+offset)
                    current_analysis.template_switch_events.append(TemplateSwitchEvent(c, i+offset, i+offset, None, None, True))
                continue
            curr_window = prev_window + c
            if len(curr_window) > n_pattern:
                curr_window = curr_window[len(curr_window) - n_pattern:]
            if curr_window in doubled_pattern:
                prev_window += c
            elif first_unit_added or len(curr_window) >= min(50, len(pattern)//2):
                if not first_unit_added:
                    first_unit_added = True
                starting_pos = []
                end_pos = []
                window_beg = prev_window[:n_pattern]
                starting_pos = self._find_all_starts(window_beg, 0) 
                if len(prev_window) > n_pattern:
                    window_end = prev_window[len(prev_window) - n_pattern:]
                    end_pos = self._find_all_starts(window_end, n_pattern - 1)
                else:
                    end_pos.append(starting_pos[0] + min(n_pattern, len(prev_window)) - 1)
                    end_pos[0] = end_pos[0] % n_pattern
                if len(starting_pos) > 2:
                    current_analysis.template_switch_event_indexes.append(i + offset)
                    current_analysis.template_switch_events.append(TemplateSwitchEvent(prev_window, i+offset-len(prev_window) + 1, i+offset, "ambiguous", "ambiguous", False))
                else:
                    current_analysis.template_switch_event_indexes.append(i + offset)
                    end_pos[0] = end_pos[0] % n_pattern
                    current_analysis.template_switch_events.append(TemplateSwitchEvent(prev_window, i+offset-len(prev_window) + 1, i+offset, starting_pos[0], end_pos[0], False))
                prev_window = c
                if prev_window not in doubled_pattern:
                    current_analysis.template_switch_event_indexes.append(i+ 1 + offset)
                    current_analysis.template_switch_events.append(TemplateSwitchEvent(prev_window, i+1+offset, i+1+offset, None, None, True))
                    prev_window = ""
            else:
                prev_window = ""

        i = len(telomer_str) - 1
        if prev_window != "" and (first_unit_added or len(prev_window) >= min(50, len(pattern)//2)):
            starting_pos = []
            end_pos = []
            window_beg = prev_window[:n_pattern]
            starting_pos = self._find_all_starts(window_beg, 0)
            if len(prev_window) > n_pattern:
                window_end = prev_window[len(prev_window) - n_pattern:]
                end_pos = self._find_all_starts(window_end, n_pattern - 1)
            else:
                end_pos.append(starting_pos[0] + min(n_pattern, len(prev_window)) - 1)
                end_pos[0] = end_pos[0] % n_pattern
            if len(starting_pos) > 2:
                current_analysis.template_switch_event_indexes.append(i + offset)
                current_analysis.template_switch_events.append(TemplateSwitchEvent(prev_window, i+offset-len(prev_window) + 1, i+offset, "ambiguous", "ambiguous", False))
            else:
                current_analysis.template_switch_event_indexes.append(i + offset)
                end_pos[0] = end_pos[0] % n_pattern
                current_analysis.template_switch_events.append(TemplateSwitchEvent(prev_window, i+offset-len(prev_window) + 1, i+offset, starting_pos[0], end_pos[0], False))
            

        return current_analysis


    def execute(self) -> None:
        #analysis_output: List[Optional[TemplateSwitchData]] = [None] * self.config.max_ends
        for telomer in self.telomers:
            if telomer and telomer.sequence:
                telomer.analysis = self._identify_template_switches(telomer.sequence)
        #return analysis_output