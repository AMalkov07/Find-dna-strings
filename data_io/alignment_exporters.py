from typing import List, Optional

from utils.data_structures import AlignmentData, Config, TelomereSequence, ImperfectAlignmentEvent

class AlignmentPrint:
    #fix figure out way to associate telomers and analysis arrays with each other
    def __init__(self, analysis: List[Optional[AlignmentData]], telomers: List[Optional[TelomereSequence]], config: Config, pattern: str):
        self.analysis = analysis
        self.telomers = telomers
        self.config = config
        self.pattern = pattern

    def _chr_end_print(self, telomer: TelomereSequence, alignment_analysis: Optional[AlignmentData], main_output_file) -> None:
        telomer_str: Optional[str] = telomer.sequence
        n_pattern = len(self.pattern)
        output = ["_______________________________"]
        output.append(telomer.chromosome_end_id)

        perfect_indexes: List[int] = alignment_analysis.perfect_alignments
        alignments: List[ImperfectAlignmentEvent] = alignment_analysis.imperfect_alignments
        all_imperfect_indexes: List[Optional[List[int]]] = [i.full_telomer_start_index for i in alignment_analysis.imperfect_alignments]
        imperfect_indexes:List[int] = [x for sublist in all_imperfect_indexes for x in sublist]

        i = j = 0
        last_end = None
        while i < len(perfect_indexes) and j < len(imperfect_indexes):
            #perfect matches
            if perfect_indexes[i] < imperfect_indexes[j]:
                if last_end and last_end < perfect_indexes[i]:
                    output.append(f"gap: {perfect_indexes[i] - last_end - 1} bp")
                last_end = perfect_indexes[i] + n_pattern
                output.append(f"{perfect_indexes[i]} (perfect match)")
                i += 1
            
            #imperfect match
            else:
                insertions_arr = alignments[j].insertion_events
                deletions_arr = alignments[j].deletion_events
                mismatches_arr = alignments[j].mismatch_events

                if last_end and last_end < imperfect_indexes[j]:
                    output.append(f"gap: {imperfect_indexes[j] - last_end - 1} bp") # check this makes sense
                last_end = imperfect_indexes[j] + n_pattern + len(insertions_arr) - len(deletions_arr)
                output.append(f"{imperfect_indexes[j]} (imperfect match), insertions: {insertions_arr}, deletions: {deletions_arr}, mismatches: {mismatches_arr}")
                j+=1

        while i < len(perfect_indexes):
            #remaining perfect matches
            if last_end and last_end < perfect_indexes[i]:
                output.append(f"gap: {perfect_indexes[i] - last_end - 1} bp")
            last_end = perfect_indexes[i] + n_pattern
            output.append(f"{perfect_indexes[i]} (perfect match)")
            i += 1

        while j < len(imperfect_indexes):
            insertions_arr = alignments[j].insertion_events
            deletions_arr = alignments[j].deletion_events
            mismatches_arr = alignments[j].mismatch_events

            if last_end and last_end < imperfect_indexes[j]:
                output.append(f"gap: {imperfect_indexes[j] - last_end - 1} bp") # check this makes sense
            last_end = imperfect_indexes[j] + n_pattern + len(insertions_arr) - len(deletions_arr)
            output.append(f"{imperfect_indexes[j]} (imperfect match), insertions: {insertions_arr}, deletions: {deletions_arr}, mismatches: {mismatches_arr}")
            j+=1
            
        print("\n".join(output), file=main_output_file) 



                

    def print_analysis(self) -> None:
        main_output_file = open(self.config.output_file, 'w')
        for i, telomer in enumerate(self.telomers):
            if telomer and telomer.sequence:
                self._chr_end_print(telomer, self.analysis[i], main_output_file)