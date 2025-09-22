from typing import List, Optional
from os.path import splitext
from _io import TextIOWrapper

from utils.data_structures import TemplateSwitchData, TemplateSwitchEvent, Config, TelomereSequence

class TemplateSwitchingPrint:
    #fix figure out way to associate telomers and analysis arrays with each other
    def __init__(self, telomers: List[Optional[TelomereSequence]], config: Config, pattern: str):
        self.telomers = telomers
        self.config = config
        self.pattern = pattern

    def _circular_distance(self, i, j, n) -> int:
        diff = abs(i - j)
        return min(diff, n-diff)  

    def _chr_end_print(self, telomer: TelomereSequence, template_switch_analysis: TemplateSwitchData, main_output_file: TextIOWrapper, variants_output_file: TextIOWrapper) -> None:
        telomer_str: str = telomer.sequence
        n_telomer_str = len(telomer_str)
        n_pattern = len(self.pattern)
        indexes = template_switch_analysis.template_switch_event_indexes

        primary_output: List[str] = ["______________________"]
        primary_output.append(telomer.chromosome_end_id)
        variants_file_output: List[str] = []

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
            strain_name: str = telomer.survivor_name
            chr_end: str = telomer.chromosome_end_id

            last_end = None
            last_last_end = None
            last_length = None

            if is_mutation:
                primary_output.append(f"pos: {telomer_start}: {telomer_chunk} mutation")
                #variants_file_output.append(f"{strain_name},{chr_end},{repeat_num},1,{reference_start},{reference_end},N/A,N/A,N/A,N/A")
                variants_file_output.append(f"{strain_name},{chr_end},,1,{telomer_start},{telomer_end},N/A,N/A,N/A,N/A")
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
                    if_small_jump = self._circular_distance(last_end, pattern_start, n_pattern) <= 5
                else:
                    if_small_jump = "N/A"
                primary_output.append(f"telomer span: {telomer_start}-{telomer_end}, length: {n_telomer_chunk}, pattern Start: {pattern_start}, pattern end: {pattern_end}")
                #variants_file_output.append(f"{strain_name},{chr_end},{repeat_num},{n_reference_chunk},{reference_start},{reference_end},{circleString_start},{circleString_end},{if_small_jump},{memory_jump_val}")
                variants_file_output.append(f"{strain_name},{chr_end},{n_telomer_chunk},{telomer_start},{telomer_end},{pattern_start},{pattern_end},{if_small_jump},{memory_jump_val}")
                last_last_end = last_end
                last_end = pattern_end
                last_length = n_telomer_chunk
            #repeat_num += 1

        final_string = "\n".join(primary_output)
        print(final_string, file=main_output_file)

        final_variant_string = "\n".join(variants_file_output)
        print(final_variant_string, file=variants_output_file)


    def print_analysis(self) -> None:
        output_file_name = self.config.output_file
        base, ext = splitext(output_file_name)
        variants_filename = f"{base}_variants.vcf"
        main_output_file = open(output_file_name, 'w')
        variants_output_file = open(variants_filename, 'w')
        for telomer in self.telomers:
            if telomer and telomer.sequence:
                if isinstance(telomer.analysis, TemplateSwitchData):
                    self._chr_end_print(telomer, telomer.analysis, main_output_file, variants_output_file)
                else:
                    raise ValueError("telomer _chr_end_print function failed to be called because the analysis type isn't TemplateSwitchData")