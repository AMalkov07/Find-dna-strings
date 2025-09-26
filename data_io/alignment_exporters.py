from typing import List, Optional, Dict, Tuple
from csv import reader
import re
from os.path import splitext

from utils.data_structures import AlignmentData, Config, TelomereSequence, ImperfectAlignmentEvent, CsvLine

class AlignmentPrint:
    def __init__(self, telomers: List[Optional[TelomereSequence]], config: Config, pattern: str):
        self.telomers = telomers
        self.config = config
        self.pattern = pattern
        self.compare_mutagenic_zone_count_match: int = 0
        self.compare_alignment_count_match: int = 0 #numer of times a mutagenic zone had the same number of matches as Ivans data
        self.exact_call_match: int = 0
        self.alignment_comparison_count: int = 0 #total number of times I compared a mutagenic arrow to Ivans data
        self.fully_matching_alignments_in_chr_end = False

    def _count_events(self, indel_arr: List[Tuple[int, str]], mode: str) -> int:
        if not indel_arr:
            return 0

        count = 1
        prev_first = indel_arr[0]

        current_event = [prev_first[1]]
        event_starts_arr = [prev_first[0]]

        for i in range(1, len(indel_arr)):
            current_first = indel_arr[i]
            condition = False
            # Only increment if not same and not consecutive
            if mode == "insertions":
                condition = (current_first[0] == prev_first[0])
            elif mode == "deletions":
                condition = (current_first[0] == prev_first[0] + 1)
            else:
                raise ValueError(f"invalid mode value: {mode}")
            
            if not (condition):
                count += 1
                current_event = [current_first[1]]
            else:
                current_event.append(current_first[1]) 
            prev_first = current_first

        return count

        
        

    def _parse_alignment_csv_line(self, line: str) -> Optional[CsvLine]:

        #pattern = re.compile(r"^IT(\d+)\s(\d+[LR])-(\d+)$")
        #pattern = re.compile(r"^KRLT(\d+)\s(\d+[LR])-(\d+)$")
        pattern = re.compile(r"^(?:IT|KRLT)(\d+)\s(\d+[LR])-(\d+)$")
        insertions = []
        deletions = []
        mismatches = []

        parts = list(reader([line]))[0]

        # Extract ID from the first field (e.g., "IT184 1L-1")
        name = parts[0]
        match = pattern.match(name)
        if match:
            survivor_id = match.group(1)
            chr_end = match.group(2)
            alignment_id = match.group(3)

            annotations = parts[1:]  # everything after the name
            for i, cell in enumerate(annotations):
                index = i + 1  # Convert to 1-based index
                cell = cell.strip()
                if not cell:
                    continue  # match
                elif cell.startswith('- '):
                    # Deletion: "- T" → deleted T from target
                    deletions.append((index, cell[2:]))
                elif cell.startswith('+ '):
                    # Insertion: "+ GTG" → inserted GTG into query
                    inserted = cell[2:]
                    for base in inserted:
                        insertions.append((index, base))
                elif ' to ' in cell:
                    # Mismatch: "T to G" → target was T, query has G
                    try:
                        from_base, to_base = map(str.strip, cell.split('to'))
                        mismatches.append((index, from_base, to_base))
                    except ValueError:
                        pass  # ignore malformed cells

            #return survivor_id, chr_end, alignment_id, insertions, deletions, mismatches
            row_data = CsvLine(
                survivor_id=int(survivor_id),
                chr_end=chr_end,
                alignment_id=int(alignment_id),
                insertions=insertions,
                deletions=deletions,
                mismatches=mismatches
            )
            return row_data
        return None

    def _sort_and_group_csv_line(self, parsed_csv_line: Dict[str, List[CsvLine]]) -> Dict[str, List[List[CsvLine]]]:
        grouped_by_chr:Dict[str, List[List[CsvLine]]] = {}
        for key in parsed_csv_line.keys():
            tmp_arr: List[CsvLine] = parsed_csv_line[key]
            tmp_arr = sorted(tmp_arr, key=lambda x: x.alignment_id)
            groups = [[tmp_arr[0]]]
            last_id = tmp_arr[0].alignment_id
            for elem in tmp_arr[1:]:
                if elem.alignment_id == last_id + 1:
                    groups[-1].append(elem)
                else:
                    groups.append([elem])
                last_id = elem.alignment_id
            grouped_by_chr[key] = groups
        return grouped_by_chr

    def _csv_parser(self, compare_file) -> Dict[str, List[List[CsvLine]]]:
        parsed_csv_line: Dict[str, List[CsvLine]] = {}
        for row in compare_file:
            line = ','.join(row)
            row_data: Optional[CsvLine] = self._parse_alignment_csv_line(line)
            if row_data:
                parsed_csv_line.setdefault(row_data.chr_end, []).append(row_data)
        grouped_by_chr: Dict[str, List[List[CsvLine]]] = self._sort_and_group_csv_line(parsed_csv_line)
        return grouped_by_chr

    def _group_new_data(self, analysis_data_grouped_by_chr: Dict[str, AlignmentData]) -> Dict[str, List[List[ImperfectAlignmentEvent]]]:
        imperfect_alignments_grouped_by_indexes: Dict[str, List[List[ImperfectAlignmentEvent]]] = {}
        for key in analysis_data_grouped_by_chr.keys():
            perfect_alignment_indexes: List[int] = sorted(analysis_data_grouped_by_chr[key].perfect_alignments)
            imperfect_alignment_indexes: List[ImperfectAlignmentEvent] = sorted(analysis_data_grouped_by_chr[key].imperfect_alignments, key=lambda event: -1 if event.full_telomer_start_index is None else event.full_telomer_start_index) # the if is None stuff is just to remove the error for the case of full_telomer_start_index being None, even thoug that should never happen
            if len(imperfect_alignment_indexes) == 0:
                imperfect_alignments_grouped_by_indexes[key] = []
                continue
            groups: List[List[ImperfectAlignmentEvent]] = []
            current_group = [imperfect_alignment_indexes[0]]
            for prev, curr in zip(imperfect_alignment_indexes, imperfect_alignment_indexes[1:]):
                # Check if any index is strictly between prev and curr
                if isinstance(prev.full_telomer_start_index, int) and isinstance(curr.full_telomer_start_index, int): 
                    if any(i > prev.full_telomer_start_index and i < curr.full_telomer_start_index for i in perfect_alignment_indexes):
                        groups.append(current_group)
                        current_group = [curr]
                    else:
                        current_group.append(curr)
                else:
                    raise ValueError("full_telomer_start_index value for some reason is None when it shouldn't be")

            groups.append(current_group)
            imperfect_alignments_grouped_by_indexes[key] = groups

        return imperfect_alignments_grouped_by_indexes
        
    #Fix: This function is prob pointless
    def _parse_analysis_data(self) -> Dict[str, List[List[ImperfectAlignmentEvent]]]:
        analysis_data_grouped_by_chr: Dict[str, AlignmentData] = {}
        for telomer in self.telomers:
            if telomer and telomer.sequence:
                chr_end = telomer.chromosome_end_id
                if isinstance(telomer.analysis, AlignmentData):
                    analysis_data_grouped_by_chr[chr_end] = telomer.analysis
                else:
                    raise ValueError("the _parse_analysis_data function failed due to telomer.analysis type not being of type telomer.analysis")
        imperfect_alignments_grouped_by_indexes = self._group_new_data(analysis_data_grouped_by_chr)
        return imperfect_alignments_grouped_by_indexes

    def _print_perfect(self, perfect_index: int):
        n_pattern = len(self.pattern)
        if self.last_end and self.last_end < perfect_index:
            self.output.append(f"gap: {perfect_index - self.last_end - 1} bp")
        self.last_end = perfect_index + n_pattern
        self.output.append(f"{perfect_index} (perfect match)")

    def _print_imperfects(self, curr_group: List[ImperfectAlignmentEvent], compare_group: Optional[List[CsvLine]], print_compare: bool) -> None:
        n_pattern = len(self.pattern)

        compare_mutagenic_zone = True
        if not print_compare or not compare_group:
            compare_mutagenic_zone = False
            self.fully_matching_alignments_in_chr_end = False
        elif len(curr_group) != len(compare_group):
            self.output.append(f"alignment number mismatch in mutagenic zone, comparison will not be made\nIvan alignment count: {len(compare_group)}, Program alignment count: {len(curr_group)}")
            compare_mutagenic_zone = False
            self.fully_matching_alignments_in_chr_end = False

        for k, val in enumerate(curr_group):
            insertions_arr: List[Tuple[int, str]] = val.insertion_events
            deletions_arr: List[Tuple[int, str]] = val.deletion_events
            mismatches_arr: List[Tuple[int, str, str]] = val.mismatch_events


            if self.last_end and self.last_end < val.full_telomer_start_index:
                self.output.append(f"gap: {val.full_telomer_start_index - self.last_end - 1} bp") # check this makes sense
            self.last_end = val.full_telomer_start_index + n_pattern + len(insertions_arr) - len(deletions_arr)
            self.output.append(f"{val.full_telomer_start_index} (imperfect match), insertions: {insertions_arr}, deletions: {deletions_arr}, mismatches: {mismatches_arr}")
                
            # code for comparing Ivans data to program self.output data
            if compare_mutagenic_zone and compare_group:
                self.alignment_comparison_count += 1
                curr_compare_group = compare_group[k]
                compare_insertions: List[Tuple[int, str]] = curr_compare_group.insertions
                compare_deletions: List[Tuple[int, str]] = curr_compare_group.deletions
                compare_mismatches: List[Tuple[int, str, str]] = curr_compare_group.mismatches

                if insertions_arr == compare_insertions and deletions_arr == compare_deletions and mismatches_arr == compare_mismatches:
                    self.output.append("alignment matches perfectly with Ivan's data")
                    self.exact_call_match += 1

                else:
                    tmp: List[str] = []
                    tmp.append(f"comparison: (Ivans Data -> self.output data): ")
                    if insertions_arr != compare_insertions:
                        tmp.append(f"Insertions: {compare_insertions} -> {insertions_arr}   ")
                    if deletions_arr != compare_deletions:
                        tmp.append(f"Deletions: {compare_deletions} -> {deletions_arr}   ")
                    if mismatches_arr != compare_mismatches:
                        tmp.append(f"Mismatches: {compare_mismatches} -> {mismatches_arr}   ")

                    self.output.append("".join(tmp))
        

    def _chr_end_print(self, telomer: TelomereSequence, alignment_analysis_data_grouped_by_chr: List[List[ImperfectAlignmentEvent]], compare_file_data_grouped_by_chr: Optional[List[List[CsvLine]]], main_output_file, stats_output_file) -> None:
        self.fully_matching_alignments_in_chr_end = True
        telomer_str: Optional[str] = telomer.sequence
        self.output: List[str] = []
        self.output.append("________________________")
        self.output.append(telomer.chromosome_end_id)
        print_compare: bool = True
        if not self.config.compare_file_path:
            print_compare = False
            self.fully_matching_alignments_in_chr_end = False
        elif not compare_file_data_grouped_by_chr or len(alignment_analysis_data_grouped_by_chr) != len(compare_file_data_grouped_by_chr):
            self.output.append(f"mutagenic zone number mismatch, comparison will not be made")
            print_compare = False
            self.fully_matching_alignments_in_chr_end = False
        else:
            self.compare_mutagenic_zone_count_match += 1
        
        self.last_end = None

        if isinstance(telomer.analysis, AlignmentData):
            perfect_indexes: List[int] = telomer.analysis.perfect_alignments

        else:
            raise ValueError("alignment _chr_end_print function failed do to telomer.analysis value not being AlignmentData")


        i = j = 0
        while i < len(perfect_indexes) and j < len(alignment_analysis_data_grouped_by_chr):
            curr_group = alignment_analysis_data_grouped_by_chr[j]
            if print_compare and compare_file_data_grouped_by_chr:
                compare_group: Optional[List[CsvLine]] = compare_file_data_grouped_by_chr[j]
            else:
                compare_group = None
            #perfect matches
            if perfect_indexes[i] < curr_group[0].full_telomer_start_index:
                self._print_perfect(perfect_indexes[i])
                i += 1
            
            #imperfect match
            else:
                self._print_imperfects(curr_group, compare_group, print_compare)
                j+=1

        while i < len(perfect_indexes):
            #remaining perfect matches
            self._print_perfect(perfect_indexes[i])
            i += 1

        while j < len(alignment_analysis_data_grouped_by_chr):
            curr_group = alignment_analysis_data_grouped_by_chr[j]
            if print_compare and compare_file_data_grouped_by_chr:
                compare_group: Optional[List[CsvLine]] = compare_file_data_grouped_by_chr[j]
            else:
                compare_group = None
            self._print_imperfects(curr_group, compare_group, print_compare)
            j+=1
            
        if self.fully_matching_alignments_in_chr_end:
            self.compare_alignment_count_match += 1
        print("\n".join(self.output), file=main_output_file) 

    def _print_stats(self, alignment_analysis_data_grouped_by_chr: Dict[str, List[List[ImperfectAlignmentEvent]]], compare_file_data_grouped_by_chr: Optional[Dict[str, List[List[CsvLine]]]], stats_output_file) -> None:
        total_chr_ends_with_mutagenic_zone = 0

        for key in alignment_analysis_data_grouped_by_chr.keys():
            if len(alignment_analysis_data_grouped_by_chr[key]) > 0 or (compare_file_data_grouped_by_chr and key in compare_file_data_grouped_by_chr and len(compare_file_data_grouped_by_chr[key]) > 0):
                total_chr_ends_with_mutagenic_zone += 1
        
        print(f"total chr ends with mutagenic zone found: {total_chr_ends_with_mutagenic_zone}", file=stats_output_file)
        if self.config.compare_file_path:
            print(f"total chr ends with matching number of mutagenic areas: {self.compare_mutagenic_zone_count_match}", file=stats_output_file)
            print(f"total chr ends with matching number of alignments: {self.compare_alignment_count_match}", file=stats_output_file)
            print(f"{self.exact_call_match} alignments matched perfectly out of {self.alignment_comparison_count} that were compared", file=stats_output_file)

        total_survivor_insertions_count = 0
        total_survivor_deletion_count = 0
        total_survivor_mismatch_count = 0
        for telomer in self.telomers:
            if telomer and isinstance(telomer.analysis, AlignmentData) and telomer.analysis.imperfect_alignments:
                chr_alignments: List[ImperfectAlignmentEvent] = telomer.analysis.imperfect_alignments
                for alignment in chr_alignments:
                    total_survivor_insertions_count += self._count_events(alignment.insertion_events, "insertions")
                    total_survivor_deletion_count += self._count_events(alignment.deletion_events, "deletions")
                    total_survivor_mismatch_count += len(alignment.mismatch_events)
        print(f"total insertion count: {total_survivor_insertions_count},  total deletions count: {total_survivor_deletion_count}, total mismatches count: {total_survivor_mismatch_count}", file=stats_output_file)


                



    def print_analysis(self) -> None:
        main_output_file = open(self.config.output_file, 'w')
        base, ext = splitext(self.config.output_file)
        stats_file_name = f"{base}_stats.txt"
        stats_output_file = open(stats_file_name, 'w')
        if self.config.compare_file_path:
            with open(self.config.compare_file_path, newline='', encoding="utf-8-sig") as csvfile:
                compare_file = reader(csvfile)
                compare_file_data_grouped_by_chr: Optional[Dict[str, List[List[CsvLine]]]] = self._csv_parser(compare_file)
        else:
            compare_file_data_grouped_by_chr = None
        alignment_analysis_data_grouped_by_chr: Dict[str, List[List[ImperfectAlignmentEvent]]] = self._parse_analysis_data()
        for telomer in self.telomers:
            if telomer and telomer.sequence:
                if isinstance(telomer.analysis, AlignmentData):
                    if compare_file_data_grouped_by_chr and telomer.chromosome_end_id in compare_file_data_grouped_by_chr.keys():
                        self._chr_end_print(telomer, alignment_analysis_data_grouped_by_chr[telomer.chromosome_end_id], compare_file_data_grouped_by_chr[telomer.chromosome_end_id], main_output_file, stats_output_file)
                    else:
                        self._chr_end_print(telomer, alignment_analysis_data_grouped_by_chr[telomer.chromosome_end_id], None, main_output_file, stats_output_file)
                else:
                    raise ValueError("the alignment _chr_end_print failed to be called because the teloemr.analysis variable type wasn't AlignmentData")
        self._print_stats(alignment_analysis_data_grouped_by_chr, compare_file_data_grouped_by_chr, stats_output_file)