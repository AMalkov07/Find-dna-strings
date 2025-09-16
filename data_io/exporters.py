from typing import List

from utils.data_structures import TemplateSwitchData, TemplateSwitchEvent, Config, TelomereSequence

class TemplateSwitchingPrint:
    #fix figure out way to associate telomers and analysis arrays with each other
    def __init__(self, analysis: List[TemplateSwitchData], telomers: TelomereSequence, config: Config, pattern: str):
        self.analysis = analysis
        self.telomers = telomers
        self.config = config
        self.pattern = pattern

    def _circular_distance(self, i, j, n) -> int:
        diff = abs(i - j)
        return min(diff, n-diff)  

    def _chr_end_print(self, telomer: TelomereSequence, template_switch_analysis: TemplateSwitchData, main_output_file) -> None:
        telomer_str: str = telomer.sequence
        n_telomer_str = len(telomer_str)
        n_pattern = len(self.pattern)
        indexes = template_switch_analysis.template_switch_event_indexes

        primary_output = [telomer.chromosome_end_id]

        for i, indx in enumerate(indexes):
            if not indx:
                continue
            template_switch_info: TemplateSwitchEvent = template_switch_analysis.template_switch_events[i]
            telomer_chunk: str = template_switch_info.telomer_chunk
            n_telomer_chunk = len(telomer_chunk)
            telomer_start: int = template_switch_info.telomer_start
            telomer_end: int = template_switch_info.telomer_end
            pattern_start: int = template_switch_info.pattern_start
            pattern_end: int = template_switch_info.pattern_end
            is_mutation: bool = template_switch_info.is_mutation

            last_end = None
            last_last_end = None
            last_length = None

            if is_mutation:
                primary_output.append(f"pos: {telomer_start}: {telomer_chunk} mutation")
                #variants_file_output.append(f"{strain_name},{chr_end},{repeat_num},1,{reference_start},{reference_end},N/A,N/A,N/A,N/A")
                last_last_end = None
                last_end = None
                last_length = None           
            else:
                if last_last_end and last_last_end != "ambiguous" and pattern_start != "ambiguous":
                   if pattern_start >= last_last_end + last_length - 1 and pattern_start <= last_last_end + last_length +1:
                       memory_jump_val = True
                   else:
                       memory_jump_val = False
                else:
                   memory_jump_val = "N/A" 
                if last_end and last_end != "ambiguous" and pattern_start != "ambiguous":
                    if_small_jump = self.circular_distance(last_end, pattern_start, n_pattern) <= 5
                else:
                    if_small_jump = "N/A"
                primary_output.append(f"telomer span: {telomer_start}-{telomer_end}, length: {n_telomer_chunk}, pattern Start: {pattern_start}, pattern end: {pattern_end}")
                #variants_file_output.append(f"{strain_name},{chr_end},{repeat_num},{n_reference_chunk},{reference_start},{reference_end},{circleString_start},{circleString_end},{if_small_jump},{memory_jump_val}")
                last_last_end = last_end
                last_end = pattern_end
                last_length = n_telomer_chunk
            #repeat_num += 1

        final_string = "\n".join(primary_output)
        print(final_string, file=main_output_file)


    def print_analysis(self) -> None:
        main_output_file = open(self.config.output_file, 'w')
        for i, telomer in enumerate(self.telomers):
            if telomer:
                self._chr_end_print(telomer, self.analysis[i], main_output_file)