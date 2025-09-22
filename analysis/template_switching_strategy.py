from utils.data_structures import TelomereSequence, Config, TemplateSwitchEvent, TemplateSwitchData
from typing import List, Optional

class TemplateSwitchingStrategy:
    def __init__(self, telomeres: List[Optional[TelomereSequence]], pattern: str, config: Config):
        self.telomers: List[Optional[TelomereSequence]] = telomeres
        self.pattern: str = pattern
        self.config = config
        
    #Fix: read through this again and break up into other functions
    def _identify_template_switches(self, telomer_str: str) -> Optional[TemplateSwitchData]:
        pattern = self.pattern
        doubled_pattern = pattern + pattern
        n_telomer = len(telomer_str)
        n_pattern = len(self.pattern)
        offset = 0
        prev_window = ""
        current_analysis = TemplateSwitchData([], [])
        for i, c in enumerate(telomer_str):
            if c == 'A' and len(prev_window) == 0:
                #self.my_dict[key].template_switching_indexes.append(i + offset)
                current_analysis.template_switch_event_indexes.append(i+offset)
                current_analysis.template_switch_events.append(TemplateSwitchEvent(c, i+offset, i+offset, None, None, True))
                continue
            curr_window = prev_window + c
            if len(curr_window) > n_pattern:
                curr_window = curr_window[len(curr_window) - n_pattern:]
            if curr_window in doubled_pattern:
                prev_window += c
            else:
                starting_pos = []
                end_pos = []
                window_beg = prev_window[:n_pattern]
                tmp_doubled_pattern= doubled_pattern
                while True:
                    start = tmp_doubled_pattern.find(window_beg) + 1
                    if start == 0:
                        break
                    starting_pos.append(start % n_pattern)
                    tmp_doubled_pattern = tmp_doubled_pattern[start + len(window_beg):]
                tmp_doubled_pattern = doubled_pattern
                if len(prev_window) > n_pattern:
                    window_end = prev_window[len(prev_window) - n_pattern:]
                    while True:
                        end = tmp_doubled_pattern.find(window_end) + 1
                        if end == 0:
                            break
                        end += n_pattern - 1
                        end_pos.append(end % n_pattern)
                        tmp_doubled_pattern = tmp_doubled_pattern[end + 1:]
                else:
                    end_pos.append(starting_pos[0] + min(n_pattern, len(prev_window)) - 1)
                    end_pos[0] %= n_pattern
                if len(starting_pos) > 2:
                    current_analysis.template_switch_event_indexes.append(i + offset)
                    current_analysis.template_switch_events.append(TemplateSwitchEvent(prev_window, i+offset-len(prev_window) + 1, i+offset, "ambiguous", "ambiguous", False))
                else:
                    current_analysis.template_switch_event_indexes.append(i + offset)
                    current_analysis.template_switch_events.append(TemplateSwitchEvent(prev_window, i+offset-len(prev_window) + 1, i+offset, starting_pos[0], end_pos[0], False))
                prev_window = c
                if prev_window not in doubled_pattern:
                    current_analysis.template_switch_event_indexes.append(i+ 1 + offset)
                    current_analysis.template_switch_events.append(TemplateSwitchEvent(prev_window, i+1+offset, i+1+offset, None, None, True))
                    prev_window = ""

        return current_analysis


    def execute(self) -> None:
        #analysis_output: List[Optional[TemplateSwitchData]] = [None] * self.config.max_ends
        for telomer in self.telomers:
            if telomer and telomer.sequence:
                telomer.analysis = self._identify_template_switches(telomer.sequence)
        #return analysis_output