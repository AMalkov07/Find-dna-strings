import csv
import re
import sys

def parse_alignment_csv_line(line: str):
    pattern = re.compile(r"^IT(\d+)\s(\d+[LR])-(\d+)$")
    #pattern = re.compile(r"^KRLT(\d+)\s(\d+[LR])-(\d+)$")
    insertions = []
    deletions = []
    mismatches = []

    parts = list(csv.reader([line]))[0]

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
        row_data = {"survivor_id": int(survivor_id),
                       "chr_end": chr_end,
                       "alignment_id": int(alignment_id),
                       "insertions": insertions,
                       "deletions": deletions,
                       "mismatches": mismatches}
        return row_data
    return None

def sort_and_group_given(grouped_by_chr_dict):
    for key in grouped_by_chr_dict.keys():
        tmp_arr = grouped_by_chr_dict[key]
        tmp_arr = sorted(tmp_arr, key=lambda x: x["alignment_id"])
        groups = [[tmp_arr[0]]]
        last_id = tmp_arr[0]["alignment_id"]
        for elem in tmp_arr[1:]:
            if elem["alignment_id"] == last_id + 1:
                groups[-1].append(elem)
            else:
                groups.append([elem])
            last_id = elem["alignment_id"]
        grouped_by_chr_dict[key] = groups
    return grouped_by_chr_dict
        



def run_csv_parser(reader):
    grouped_by_chr_dict = {}
    for row in reader:
        line = ','.join(row)
        row_data = parse_alignment_csv_line(line)
        if row_data:
            key = row_data["chr_end"]
            grouped_by_chr_dict.setdefault(key, []).append(row_data)
    grouped_by_chr_dict = sort_and_group_given(grouped_by_chr_dict)
    return grouped_by_chr_dict

def group_new_data(new_data_grouped_by_chr_dict):
    for key in new_data_grouped_by_chr_dict.keys():
        indexes = sorted(new_data_grouped_by_chr_dict[key].indexes)
        extra_indexes = sorted(new_data_grouped_by_chr_dict[key].extra_alignment_indexes)
        if len(extra_indexes) == 0:
            continue
        groups = []
        current_group = [extra_indexes[0]]
        for prev, curr in zip(extra_indexes, extra_indexes[1:]):
            # Check if any index is strictly between prev and curr
            if any(i > prev and i < curr for i in indexes):
                groups.append(current_group)
                current_group = [curr]
            else:
                current_group.append(curr)

        groups.append(current_group)
        new_data_grouped_by_chr_dict[key].extra_alignment_indexes = groups
    return new_data_grouped_by_chr_dict
    

def parse_new_data(new_data):
    new_data_grouped_by_chr_dict = {}
    pattern = re.compile(r"^([^_]+)_(\d+[LR])")
    for elem in new_data:
        name = elem[0]
        match = pattern.match(name)
        if not match:
            print("no match was found")
        strain_num, chr_end = match.groups()
        #new_data_grouped_by_chr_dict[match.group(1)] = elem[1]
        new_data_grouped_by_chr_dict[chr_end] = elem[1]
    new_data_grouped_by_chr_dict = group_new_data(new_data_grouped_by_chr_dict)
    return new_data_grouped_by_chr_dict, strain_num

def count_tuples(arr, mode):
    if not arr:
        return 0
    
    count = 1  # start by counting the first tuple
    prev_first = arr[0][0]

    for i in range(1, len(arr)):
        current_first = arr[i][0]
        condition = False
        # Only increment if not same and not consecutive
        if mode == "insertions":
            condition = (current_first == prev_first)
        elif mode == "deletions":
            condition = (current_first == prev_first + 1)
            
        if not (condition):
            count += 1
        prev_first = current_first

    return count

def count_events(arr, mode):
    if not arr:
        return 0, 0, 0

    total_events = []
    
    count = 1  # start by counting the first tuple
    prev_first = arr[0]

    current_event = [prev_first[1]]
    ref_start = [prev_first[0]]

    for i in range(1, len(arr)):
        current_first = arr[i]
        condition = False
        # Only increment if not same and not consecutive
        if mode == "insertions":
            condition = (current_first[0] == prev_first[0])
        elif mode == "deletions":
            condition = (current_first[0] == prev_first[0] + 1)
            
        if not (condition):
            count += 1
            total_events.append("".join(current_event))
            current_event = [current_first[1]]
            ref_start.append(current_first[0])
        else:
            current_event.append(current_first[1]) 
        prev_first = current_first
    total_events.append("".join(current_event))

    return count, total_events, ref_start

#note: self is just cuz I'm too lazy to change the variable name, it's not actually an object function
def custom_object_print(self, given_data_grouped_by_chr_dict, print_comparison, strain_num, chr_end, variants_filename, referenceString):
    total_insertions = 0
    total_deletions = 0
    total_mismatches = 0
    total_given_insertions = 0
    total_given_deletions = 0
    total_given_mismatches = 0
    repeat_num = 0
    alignment_mismatches_comparison = None
    if not self.extra_alignment_indexes:
        self.extra_alignment_indexes = []
    all_imperfect_alignments = []
    output = []
    total_alignments = 0
    if_matching_number_of_alignment = True
    total_number_alignments_compared_in_chr = 0
    total_number_perfect_alignment_matches_in_chr = 0
    if self.extra_alignment_indexes:
        for group in self.extra_alignment_indexes:
            total_alignments += len(group)
        total_matches = len(self.indexes) + total_alignments
    else:
        total_matches = len(self.indexes)
    output.append(f"Total matches: {total_matches}")

    # Two-pointer merge
    i = j = 0
    last_end = None
    while i < len(self.indexes) and j < len(self.extra_alignment_indexes):
        repeat_num += 1
        curr_group = self.extra_alignment_indexes[j]
        # perfect match
        if self.indexes[i] <= curr_group[0]:
            if last_end and last_end < self.indexes[i]:
                output.append(f"gap: {self.indexes[i] - last_end - 1} bp")
            last_end = self.indexes[i] + self.n_subStr
            output.append(f"{self.indexes[i]} (perfect match)")
            i += 1
        # imperfect match
        else:
            should_output_extra_alignment_error_message = True
            should_increment_repeat_num = 1
            for k in range(len(curr_group)):
                val = curr_group[k]
                insertions_arr = self.extra_alignment_insertions_and_deletions[val][0]
                deletions_arr = self.extra_alignment_insertions_and_deletions[val][1]
                mismatches_arr = self.extra_alignment_insertions_and_deletions[val][2]
                #total_insertions += len(insertions_arr)
                #total_deletions += len(deletions_arr)
                if len(insertions_arr) > 0:
                    curr_insertions, insertions_events, ref_start = count_events(insertions_arr, "insertions")
                    total_insertions += curr_insertions
                    for ii, ins in enumerate(insertions_events):
                        curr_ref_start = ref_start[ii]
                        ref_subStr = referenceString[curr_ref_start-1:curr_ref_start+1]
                        print(f"{strain_num}, {chr_end}, {repeat_num}, insertion, {len(ins)}, {ref_subStr}, {ins}, {all(c in "TG" for c in ins)}", file=variants_filename)
                if len(deletions_arr) > 0:
                    curr_deletions, deletions_events, ref_start = count_events(deletions_arr, "deletions")
                    total_deletions += curr_deletions
                    for ii, dele in enumerate(deletions_events):
                        n_dele = len(dele)
                        curr_ref_start = ref_start[ii]
                        ref_subStr = referenceString[curr_ref_start-2:curr_ref_start+n_dele]
                        print(f"{strain_num}, {chr_end}, {repeat_num}, deletion, {len(dele)}, {ref_subStr}, {dele}, {all(c in "TG" for c in dele)}", file=variants_filename)
                total_mismatches += len(mismatches_arr)
                for mis in mismatches_arr:
                    print(f"{strain_num}, {chr_end}, {repeat_num}, single base, {mis[1]} -> {mis[2]}, {mis[2] in "TG"}", file=variants_filename)

                if should_increment_repeat_num < len(curr_group):
                    should_increment_repeat_num += 1
                    repeat_num += 1
                    
                if last_end and last_end < val:
                    output.append(f"gap: {val - last_end - 1} bp") # check this makes sense
                last_end = val + self.n_subStr + len(insertions_arr) - len(deletions_arr)

                if print_comparison and should_output_extra_alignment_error_message:
                    if len(curr_group) != len(given_data_grouped_by_chr_dict[j]):
                        output.append("alignments number mismatch in current mutagenic area. No further comparison will be made for chr end")
                        if_matching_number_of_alignment = False
                        output.append(f"Ivans number of alignments: {len(given_data_grouped_by_chr_dict[j])}, output number of alignmnets: {len(curr_group)}")
                        #alignment_mismatches_comparison = (len(given_data_grouped_by_chr_dict[j]), len(curr_group))
                        alignment_mismatches_comparison = f"{len(given_data_grouped_by_chr_dict[j])} -> {len(curr_group)}"
                        should_output_extra_alignment_error_message = False

                output.append(
                    f"{val} (imperfect match), insertions: {insertions_arr}, deletions: {deletions_arr}, mismatches: {mismatches_arr}")
                if print_comparison and if_matching_number_of_alignment:
                    total_number_alignments_compared_in_chr += 1
                    if k >= len(given_data_grouped_by_chr_dict[j]):
                        output.append("ERROR More alignments in the mutagenic area then in Ivan's data")
                        continue
                    given_insertions_arr = given_data_grouped_by_chr_dict[j][k]["insertions"]
                    given_deletions_arr = given_data_grouped_by_chr_dict[j][k]["deletions"]
                    given_mismatches_arr = given_data_grouped_by_chr_dict[j][k]["mismatches"]

                    curr_given_insertions, x, y = count_events(given_insertions_arr, "insertions")
                    total_given_insertions += curr_given_insertions
                    curr_given_deletions, x, y = count_events(given_deletions_arr, "deletions")
                    total_given_deletions += curr_given_deletions
                    total_given_mismatches += len(given_mismatches_arr)
                    

                    if given_insertions_arr == insertions_arr and given_deletions_arr == deletions_arr and given_mismatches_arr == mismatches_arr:
                        output.append("alignment matches perfectly with Ivan's data ")
                        total_number_perfect_alignment_matches_in_chr +=  1
                    
                    else:
                        tmp = []
                        tmp.append(f"comparison: (Ivans Data -> output data): ")
                        if insertions_arr != given_insertions_arr:
                            tmp.append(f"Insertions: {given_insertions_arr} -> {insertions_arr}   ")
                        if deletions_arr != given_deletions_arr:
                            tmp.append(f"deletions: {given_deletions_arr} -> {deletions_arr}   ")
                        if mismatches_arr != given_mismatches_arr:
                            tmp.append(f"mismatches: {given_mismatches_arr} -> {mismatches_arr}   ")
                        ivans_events = count_tuples(given_insertions_arr, "insertions") + count_tuples(given_deletions_arr, "deletions") + count_tuples(given_mismatches_arr, "mismatches")
                        output_events = count_tuples(insertions_arr, "insertions") + count_tuples(deletions_arr, "deletions") + count_tuples(mismatches_arr, "mismatches")
                        tmp.append(f"number of events: {ivans_events} -> {output_events}")


                        output.append("".join(tmp))
                        all_imperfect_alignments.append("".join(tmp))

                # f"{self.extra_alignment_indexes[j]} (imperfect match)")
            j += 1

    # Remaining perfect matches
    while i < len(self.indexes):
        if last_end and last_end < self.indexes[i]:
            output.append(f"gap: {self.indexes[i] - last_end - 1} bp")
        last_end = self.indexes[i] + self.n_subStr
        output.append(f"{self.indexes[i]} (perfect match)")
        i += 1

    # Remaining imperfect matches
    while j < len(self.extra_alignment_indexes):
        curr_group = self.extra_alignment_indexes[j]
        should_output_extra_alignment_error_message = True
        should_increment_repeat_num = 1
        for k in range(len(curr_group)):
            val = curr_group[k]
            insertions_arr = self.extra_alignment_insertions_and_deletions[val][0]
            deletions_arr = self.extra_alignment_insertions_and_deletions[val][1]
            mismatches_arr = self.extra_alignment_insertions_and_deletions[val][2]
            #total_insertions += len(insertions_arr)
            #total_deletions += len(deletions_arr)
            #total_mismatches += len(mismatches_arr)
            if len(insertions_arr) > 0:
                curr_insertions, insertions_events, ref_start = count_events(insertions_arr, "insertions")
                total_insertions += curr_insertions
                for ii, ins in enumerate(insertions_events):
                    curr_ref_start = ref_start[ii]
                    ref_subStr = referenceString[curr_ref_start-1:curr_ref_start+1]
                    print(f"{strain_num}, {chr_end}, {repeat_num}, insertion, {len(ins)}, {ref_subStr}, {ins}, {all(c in "TG" for c in ins)}", file=variants_filename)
            if len(deletions_arr) > 0:
                curr_deletions, deletions_events, ref_start = count_events(deletions_arr, "deletions")
                total_deletions += curr_deletions
                for ii, dele in enumerate(deletions_events):
                    n_dele = len(dele)
                    curr_ref_start = ref_start[ii]
                    ref_subStr = referenceString[curr_ref_start-2:curr_ref_start+n_dele]
                    print(f"{strain_num}, {chr_end}, {repeat_num}, deletion, {len(dele)}, {ref_subStr}, {dele}, {all(c in "TG" for c in dele)}", file=variants_filename)
            total_mismatches += len(mismatches_arr)
            for mis in mismatches_arr:
                print(f"{strain_num}, {chr_end}, {repeat_num}, single base, {mis[1]} -> {mis[2]}, {mis[2] in "TG"}", file=variants_filename)

            if should_increment_repeat_num < len(curr_group):
                should_increment_repeat_num += 1
                repeat_num += 1

            if last_end and last_end < val:
                output.append(f"gap: {val - last_end - 1} bp") # check this makes sense
            last_end = val + self.n_subStr + len(insertions_arr) - len(deletions_arr)

            if print_comparison and should_output_extra_alignment_error_message:
                if len(curr_group) != len(given_data_grouped_by_chr_dict[j]):
                    output.append("alignments number mismatch in current mutagenic area. No further comparison will be made for chr end")
                    if_matching_number_of_alignment = False
                    output.append(f"Ivans number of alignments: {len(given_data_grouped_by_chr_dict[j])}, output number of alignmnets: {len(curr_group)}")
                    alignment_mismatches_comparison = f"{len(given_data_grouped_by_chr_dict[j])} -> {len(curr_group)}"
                should_output_extra_alignment_error_message = False
                

            output.append(
                # f"{val} (imperfect match), insertions: {[x + int(val) for x in self.extra_alignment_insertions_and_deletions[val][0]]}, deletions: {[x + int(val) for x in self.extra_alignment_insertions_and_deletions[val][1]]}, mismatches: {[x + int(val) for x in self.extra_alignment_insertions_and_deletions[val][2]]}")
                f"{val} (imperfect match), insertions: {insertions_arr}, deletions: {deletions_arr}, mismatches: {mismatches_arr}")
            if print_comparison and if_matching_number_of_alignment:
                total_number_alignments_compared_in_chr += 1
                if k >= len(given_data_grouped_by_chr_dict[j]):
                    output.append("More alignments in the mutagenic area then in Ivan's data")
                    continue
                given_insertions_arr = given_data_grouped_by_chr_dict[j][k]["insertions"]
                given_deletions_arr = given_data_grouped_by_chr_dict[j][k]["deletions"]
                given_mismatches_arr = given_data_grouped_by_chr_dict[j][k]["mismatches"]

                curr_given_insertions, x, y = count_events(given_insertions_arr, "insertions")
                total_given_insertions += curr_given_insertions
                curr_given_deletions, x, y = count_events(given_deletions_arr, "deletions")
                total_given_deletions += curr_given_deletions
                total_given_mismatches += len(given_mismatches_arr)

                if given_insertions_arr == insertions_arr and given_deletions_arr == deletions_arr and given_mismatches_arr == mismatches_arr:
                    output.append("alignment matches perfectly with Ivan's data ")
                    total_number_perfect_alignment_matches_in_chr +=  1
                
                else:
                    tmp = []
                    tmp.append(f"comparison: (Ivans Data -> output data): ")
                    if insertions_arr != given_insertions_arr:
                        tmp.append(f"Insertions: {given_insertions_arr} -> {insertions_arr}   ")
                    if deletions_arr != given_deletions_arr:
                        tmp.append(f"deletions: {given_deletions_arr} -> {deletions_arr}   ")
                    if mismatches_arr != given_mismatches_arr:
                        tmp.append(f"mismatches: {given_mismatches_arr} -> {mismatches_arr}")
                    ivans_events = count_tuples(given_insertions_arr, "insertions") + count_tuples(given_deletions_arr, "deletions") + count_tuples(given_mismatches_arr, "mismatches")
                    output_events = count_tuples(insertions_arr, "insertions") + count_tuples(deletions_arr, "deletions") + count_tuples(mismatches_arr, "mismatches")
                    tmp.append(f"number of events: {ivans_events} -> {output_events}")

                    output.append("".join(tmp))
                    all_imperfect_alignments.append("".join(tmp))
            # f"{self.extra_alignment_indexes[j]} (imperfect match)")
        j += 1

    final_string = "\n".join(output)
    print(final_string)

    return if_matching_number_of_alignment, total_number_alignments_compared_in_chr, total_number_perfect_alignment_matches_in_chr, ("\n".join(all_imperfect_alignments)), total_insertions, total_deletions, total_mismatches, total_given_insertions, total_given_deletions, total_given_mismatches, alignment_mismatches_comparison
    #return total_insertions, total_deletions, total_mismatches


def print_differences(new_data_grouped_by_chr_dict, given_data_grouped_by_chr_dict, stats_filename, variants_filename, strain_num, referenceString):

    print("strain name, chr end, repeat number, mutation type, mutation length, mutated area (w/ before and after bp), mutation bases, only TG mutations", file=variants_filename)

    total_matching_number_of_gaps = 0
    total_matching_number_of_alignments = 0 # the number of times that the number of alignments is the same inside of a mutagenic area
    total_number_alignments_compared = 0
    total_number_perfect_alignment_matches = 0
    total_chr_ends_compared = 0
    total_no_mutagenic_zones = 0
    total_chr_ends_with_mutagenic_zone = 0
    total_exact_match_for_entire_chr_end = 0

    mismatch_gap_number_indexes = []
    mismatch_alignment_number_indexes = []
    all_imperfect_alignments = []

    total_insertions = 0
    total_deletions = 0
    total_mismatches = 0
    total_given_insertions = 0
    total_given_deletions = 0
    total_given_mismatches = 0
    absolute_insertions = 0
    absolute_deletions = 0
    absolute_mismatches = 0

    chr_end_insertiong_deletion_count = []

    for key in new_data_grouped_by_chr_dict.keys():
        chr_end_insertion = 0
        chr_end_deletion = 0
        chr_end_mismatch = 0
        for group in new_data_grouped_by_chr_dict[key].extra_alignment_indexes:
            for elem in group:
                curr_insertions, _, _ = count_events(new_data_grouped_by_chr_dict[key].extra_alignment_insertions_and_deletions[elem][0], "insertions")
                absolute_insertions += curr_insertions
                chr_end_insertion += curr_insertions
                curr_deletions, _, _ = count_events(new_data_grouped_by_chr_dict[key].extra_alignment_insertions_and_deletions[elem][1], "deletions")
                absolute_deletions += curr_deletions
                chr_end_deletion += curr_deletions
                absolute_mismatches += len(new_data_grouped_by_chr_dict[key].extra_alignment_insertions_and_deletions[elem][2])
                chr_end_mismatch += len(new_data_grouped_by_chr_dict[key].extra_alignment_insertions_and_deletions[elem][2])
        chr_end_insertiong_deletion_count.append(f"{key}: insertions: {chr_end_insertion}, deletions: {chr_end_deletion}, mismatches: {chr_end_mismatch}")



    for key in given_data_grouped_by_chr_dict.keys():
        for i in ((given_data_grouped_by_chr_dict[key])):
            for j in ((i)):
                #print(j)
                if j:
                    curr_given_insertions, x, y = count_events(j['insertions'], "insertions")
                    total_given_insertions += curr_given_insertions
                    curr_given_deletions, x, y = count_events(j['deletions'], "deletions")
                    total_given_deletions += curr_given_deletions
                    total_given_mismatches += len(j['mismatches'])
            


    for key in new_data_grouped_by_chr_dict.keys():
        all_alignment_mismatch_comparisons = []
        total_chr_ends_compared += 1
        if len(new_data_grouped_by_chr_dict[key].extra_alignment_indexes) > 0 or (key in given_data_grouped_by_chr_dict and len(given_data_grouped_by_chr_dict[key]) > 0):
            total_chr_ends_with_mutagenic_zone += 1
        print_comparison = True
        print("______________________________________________________________")
        print(f"{key}:")
        if key not in given_data_grouped_by_chr_dict:
            if len(new_data_grouped_by_chr_dict[key].extra_alignment_indexes) > 0:
                total_alignments = 0
                for elem in new_data_grouped_by_chr_dict[key].extra_alignment_indexes:
                    total_alignments += len(elem)
                print(f"number of mutagenic zones mismatch. No alignments were found in Ivans Data. Comparison will not be made\nIvans data: 0 alignments, Output data: {total_alignments} alignments")
                mismatch_gap_number_indexes.append(key)
            #else:
                #total_matching_number_of_gaps += 1
                #total_no_mutagenic_zones += 1
            print_comparison = False
        if print_comparison:
            alignments_groups_found = new_data_grouped_by_chr_dict[key].extra_alignment_indexes
            alignment_groups_given = given_data_grouped_by_chr_dict[key]
            n_alignments_groups_found = len(alignments_groups_found)
            n_alignments_groups_given = len(alignment_groups_given)
            if n_alignments_groups_found != n_alignments_groups_given:
                print(f"number of mutagenic zones mismatch. Comparison will not be made")
                mismatch_gap_number_indexes.append(key)
                total_alignments_found = 0
                total_alignments_given = 0
                for elem in new_data_grouped_by_chr_dict[key].extra_alignment_indexes:
                    total_alignments_found += len(elem)
                for elem in alignment_groups_given:
                    total_alignments_given += len(elem)
                print(f"Ivan mutagenic zones: {n_alignments_groups_given}, output mutagenic zones: {n_alignments_groups_found}, Ivan total alignments: {total_alignments_given}, Output total alignments: {total_alignments_found}")
                print_comparison = False
            else:
                total_matching_number_of_gaps += 1
        if print_comparison:
            if_matching_number_of_alignments, total_number_alignments_compared_in_chr, total_number_perfect_alignment_alignmen_in_chr, tmp_imperfect_alignments, n_insertions, n_deletions, n_mismatches, n_given_insertions, n_given_deletions, n_given_mismatches, alignment_mismatches_comparison = custom_object_print(new_data_grouped_by_chr_dict[key], given_data_grouped_by_chr_dict[key], print_comparison, strain_num, key, variants_filename, referenceString)
            total_insertions += n_insertions
            total_deletions += n_deletions
            total_mismatches += n_mismatches
            #total_given_insertions += n_given_insertions
            #total_given_deletions += n_given_deletions
            #total_given_mismatches += n_given_mismatches
            if if_matching_number_of_alignments and len(tmp_imperfect_alignments) == 0:
                total_exact_match_for_entire_chr_end += 1
            else:
                all_imperfect_alignments.append(tmp_imperfect_alignments)
            if if_matching_number_of_alignments:
                total_matching_number_of_alignments += 1
            else:
                mismatch_alignment_number_indexes.append(str(key) + " " + alignment_mismatches_comparison)
            total_number_alignments_compared += total_number_alignments_compared_in_chr
            total_number_perfect_alignment_matches += total_number_perfect_alignment_alignmen_in_chr
        else:
            #insertions, deletions, mismatches = custom_object_print(new_data_grouped_by_chr_dict[key], [], print_comparison)
            if_matching_number_of_alignments, total_number_alignments_compared_in_chr, total_number_perfect_alignment_alignmen_in_chr, tmp_imperfect_alignments, n_insertions, n_deletions, n_mismatches, n_given_insertions, n_given_deletions, n_given_mismatches, alignment_mismatches_comparison = custom_object_print(new_data_grouped_by_chr_dict[key], [], print_comparison, strain_num, key, variants_filename, referenceString)
            #total_insertions += n_insertions
            #total_deletions += n_deletions
            #total_mismatches += n_mismatches
            total_given_insertions += n_given_insertions
            total_given_deletions += n_given_deletions
            total_given_mismatches += n_given_mismatches
        print("\n")

    # check if key is in Ivans data but not in my output data:
    for key in given_data_grouped_by_chr_dict:
        if key not in new_data_grouped_by_chr_dict.keys():
            total_chr_ends_compared += 1
            total_chr_ends_with_mutagenic_zone += 1
            print(f"{key}:")
            print(f"chr end present in Ivan's output but not present in program output")
            print("\n")

    print(f"total chr ends with mutagenic zone found: {total_chr_ends_with_mutagenic_zone}", file=stats_filename)
    print(f"total chr ends with matching number of mutagenic areas: {total_matching_number_of_gaps}", file=stats_filename)
    print(f"total chr ends with matching number of alignments: {total_matching_number_of_alignments}", file=stats_filename)
    print(f"total chr ends with exact same call as Ivans data: {total_exact_match_for_entire_chr_end}", file=stats_filename)
    print(f"{total_number_perfect_alignment_matches} alignments matched perfectly out of {total_number_alignments_compared} that were compared", file=stats_filename)
    #print(f"total insertion count: {total_insertions},  total deletions count: {total_deletions}, total mismatches count: {total_mismatches}", file=stats_filename)
    print(f"total insertion count: {absolute_insertions},  total deletions count: {absolute_deletions}, total mismatches count: {absolute_mismatches}", file=stats_filename)
    print(f"total Ivan insertion count: {total_given_insertions},  total Ivan deletions count: {total_given_deletions}, total Ivan, mismatches count: {total_given_mismatches}", file=stats_filename)
    print("\n".join(chr_end_insertiong_deletion_count), file=stats_filename)
    print("\nchr ends with mutagenic zone number mismatch:", file=stats_filename)
    print("\n".join(mismatch_gap_number_indexes), file=stats_filename)
    print("\nchr ends with alignment number mismatch:", file=stats_filename)
    print("\n".join(mismatch_alignment_number_indexes), file=stats_filename)
    print("\nall Alignment mismatches:", file=stats_filename)        
    print("\n\n".join(all_imperfect_alignments), file=stats_filename)
        
            

        

    
    
        

def compare_outputs(previous_data_csvfile, new_data, stats_filename, variants_filename, referenceString):
    with open(previous_data_csvfile, newline='', encoding="utf-8-sig") as csvfile:
        reader = csv.reader(csvfile)
        given_data_grouped_by_chr_dict = run_csv_parser(reader) # structure of this dict: {"1L": [[dict1], [d2], [d3, d4, d5 d6]], "2L": ...}
        new_data_grouped_by_chr_dict, strain_num = parse_new_data(new_data) # structure of this dict: {"1L": self_search_obj(1L), "2L": self_search_obj(2L), ...}
        print_differences(new_data_grouped_by_chr_dict, given_data_grouped_by_chr_dict, stats_filename, variants_filename, strain_num, referenceString)
        for key in new_data_grouped_by_chr_dict.keys():
            if key not in given_data_grouped_by_chr_dict:
                print(f"{key} key missing from Ivan's data")
                continue
            alignments_groups_found = new_data_grouped_by_chr_dict[key].extra_alignment_indexes
            alignment_groups_given = given_data_grouped_by_chr_dict[key]
            n_alignments_groups_found = len(new_data_grouped_by_chr_dict[key].extra_alignment_indexes)
            n_alignments_groups_given = len(given_data_grouped_by_chr_dict[key])
            return
            

#with open("181_Ivan_data.csv", newline='', encoding="utf-8-sig") as csvfile:
    #reader = csv.reader(csvfile)
    #grouped_by_chr_dict = run_csv_parser(reader)
    #for key in grouped_by_chr_dict.keys():
        #for group in grouped_by_chr_dict[key]:
            #for row_data in group:
                #survivor_id = row_data["survivor_id"]
                #chr_end = row_data["chr_end"]
                #alignment_id = row_data["alignment_id"]
                #insertions = row_data["insertions"]
                #deletions = row_data["deletions"]
                #mismatches = row_data["mismatches"]
                
                #print(f"survivor_id: {survivor_id}")
                #print(f"chr end: {chr_end}")
                #print(f"Alignment ID: {alignment_id}")
                #print("Insertions:", insertions)
                #print("Deletions:", deletions)
                #print("Mismatches:", mismatches)
                #print("-" * 40)
