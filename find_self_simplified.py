from collections import deque, defaultdict
import sys
from Bio import Align
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

fig, ax = plt.subplots(figsize=(16,4))

global_figure_name = "IT156_graph.png"
global_counter = -1
global_max_x = 0
global_min_x = 0

class filter_options:
    def __init__(self, only_non_overlapping = True, min_substring_length = 5, min_matches = 1, is_best_only = False, required_substring_length = -1):
        self.only_non_overlapping = only_non_overlapping
        self.min_substring_length = min_substring_length
        self.min_matches = min_matches
        self.is_best_only = is_best_only
        self.required_substring_length = required_substring_length

class self_search_type:
    def __init__(self, n_subStr, indexes, min_gap, is_overlap):
        self.n_subStr = n_subStr
        self.indexes = indexes
        self.min_gap = min_gap
        self.is_overlap = is_overlap
        self.n_extra_alignment = 0
        self.extra_alignment_indexes = []
        self.extra_alignment_insertions_and_deletions = {}
    
    def average_distance(self):
        arr = self.indexes
        tmp = 0
        n = len(arr)
        for i in range(1, n):
            tmp += (arr[i] - arr[i-1])
        return tmp//n

    
    def __repr__(self):
        n = self.n_subStr
        #avg_dist = self.average_distance()
        return f"substring length={self.n_subStr}, number of matches={len(self.indexes)}, min_gap={self.min_gap}, overlap={self.is_overlap}, number of extra alignments={self.n_extra_alignment}\nstarting indexes: {self.indexes}\nextra alignemtn indexes: {self.extra_alignment_indexes}\n"


aligner = Align.PairwiseAligner()

aligner.mode = 'local'  # Local (Smith-Waterman) alignment; use 'global' for global (Needleman-Wunsch)
aligner.match_score = 1  # Score for a match
aligner.mismatch_score = -(1)  # Penalty for a mismatch
aligner.open_gap_score = -(1)  # Penalty for opening a gap
aligner.extend_gap_score = -(.5)  # Penalty for extending a gap
#aligner.query_left_open_gap_score = 0
#aligner.query_right_open_gap_score = 0
#aligner.target_left_open_gap_score = 0
#aligner.target_right_open_gap_score = 0

def alignment(input_s, key, my_dict_entry):
    alignments = my_dict_entry.indexes
    key_n = len(key)

    last_val = 0
    str_dict = {}
    queue = deque()
    for alignment in alignments:
        tmp = input_s[last_val:alignment]
        str_dict[tmp] = [last_val, alignment]
        queue.append(tmp)
        last_val = alignment+key_n

    good_alignment_arr = []
    good_alignment_starting_pos = []
    dict_insertions_and_deletions = {}
    arr = []
    counter = 0
    while queue:
        str = queue.popleft()
        if str != "":
            alignments = aligner.align(str, key)
            # think about why this if else statement is necesary<<<<<
            # try except block is necessary cause of overflow if there are too many optimal alignments
            try: 
                n_alignments = len(alignments)
            except OverflowError:
                n_alignments = 1

            if n_alignments > 0:
                score = alignments[0].score
            else:
                continue
            # below if statement determines what qualifies as a good alignment score
            if score >= key_n - key_n/10:
                counter += 1 
                good_alignment_arr.append(alignments[0])
                
                alignment_start = alignments[0].coordinates[0][0]
                good_alignment_starting_pos.append(str_dict[str][0]+alignment_start)

                new_entry_start = str_dict[str][0]
                new_entry_end = new_entry_start + alignment_start
                new_str = input_s[new_entry_start:new_entry_end+1]
                str_dict[new_str] = [new_entry_start, new_entry_end]
                queue.append(new_str)

                #new_entry_start = str_dict[str][0] + alignments[0].coordinates[0][0] + key_n + 1
                new_entry_start = str_dict[str][0] + alignments[0].coordinates[0][-1] + 1
                new_entry_end = str_dict[str][1]
                new_str = input_s[new_entry_start:new_entry_end]
                str_dict[new_str] = [new_entry_start, new_entry_end]
                queue.append(new_str)

                coordinates = alignments[0].coordinates
                insertions_indexes = []
                deletions_indexes = []

                for i in range(coordinates[1][0]):
                    insertions_indexes.append(i+1)

                # Iterate over each alignment block in coordinates
                for i in range(1, len(coordinates[0])-2, 2):
                    # Extract the end of the current block and the start of the next
                    seq1_end, seq1_next = coordinates[0][i], coordinates[0][i+1]
                    seq2_end, seq2_next = coordinates[1][i], coordinates[1][i+1]

                    # Count gaps between alignment blocks
                    if seq1_end != seq1_next:  # Deletion in seq2
                        for i in range(seq1_next-seq1_end):
                           deletions_indexes.append(seq2_end + i)
                    if seq2_end != seq2_next:  # Insertion in seq1
                        for i in range(seq2_next-seq2_end):
                           insertions_indexes.append(seq2_end + i)

                dict_insertions_and_deletions[str_dict[str][0]+alignment_start] = [insertions_indexes, deletions_indexes]

            arr.append(score)
    return(good_alignment_starting_pos, arr, dict_insertions_and_deletions)

def graph_setup():
    #ax.text(0, 8, "IT148", ha='center', va='center', fontsize=12, color='black')
    #ax.text(0, 7, "chr1", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 6.5, "chr2", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 6, "chr3", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 5.5, "chr4", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 5, "chr5", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 4.5, "chr6", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 4, "chr7", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 3.5, "chr8", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 3, "chr9", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 2.5, "chr10", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 2, "chr11", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 1.5, "chr12", ha='center', va='center', fontsize=10, color='black')
    #ax.text(0, 1, "chr13", ha='center', va='center', fontsize=10, color='black')

    #ax.text(0, 9, "IT148", ha='center', va='center', fontsize=12, color='black')
    ax.text(0, 9, "IT156", ha='center', va='center', fontsize=12, color='black')
    ax.text(0, 8.5, "chr1", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 8, "chr2", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 7.5, "chr3", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 7, "chr4", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 6.5, "chr5", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 6, "chr6", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 5.5, "chr7", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 5, "chr8", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 4.5, "chr9", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 4, "chr10", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 3.5, "chr11", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 3, "chr12", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 2.5, "chr13", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 2, "chr14", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 1.5, "chr15", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 1, "chr16", ha='center', va='center', fontsize=10, color='black')

    #ax.set_xlabel("X-axis")
    #ax.set_ylim(.5, 8.5)
    ax.set_ylim(.5, 9.5)

    # Remove x/y-axis and add legend
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)

def graph_output(my_dict):
    offset = 300 
    counter = 0
    global global_counter
    global global_max_x
    global global_min_x
    #y_index = (14 - (global_counter//2)) / 2
    y_index = (17 - (global_counter//2)) / 2
    key = list(my_dict.keys())[0]
    arrow_distance_orig = my_dict[key].n_subStr-1
    offset_offset = arrow_distance_orig//2
    set1 = my_dict[key].indexes
    set2 = my_dict[key].extra_alignment_indexes
    set2.sort()
    insertions_and_deletions = my_dict[key].extra_alignment_insertions_and_deletions
    all_points = []
    set1_ptr = 0
    set2_ptr = 0
    set1_n = len(set1)
    set2_n = len(set2)
    while set1_ptr < set1_n and set2_ptr < set2_n:
        set1_tmp = set1[set1_ptr]
        set2_tmp = set2[set2_ptr]
        if set1_tmp < set2_tmp:
            #all_points.append(set1_tmp + offset)
            all_points.append(set1_tmp + offset + (offset_offset*counter))
            set1_ptr += 1
        else:
            #all_points.append(set2_tmp + offset)
            all_points.append(set2_tmp + offset + (offset_offset*counter))
            set2_ptr += 1
        counter += 1
    if set1_ptr == set1_n:
        for i in range(set2_ptr, set2_n):
            #all_points.append(set2[i] + offset)
            all_points.append(set2[i] + offset + (offset_offset* counter))
            counter += 1
    else:
        for i in range(set1_ptr, set1_n, 1):
            #all_points.append(set1[i] + offset)
            all_points.append(set1[i] + offset + (offset_offset*counter))
            counter += 1
    
    #arrow_style = "->"
    sign = global_counter % 2
    if sign == 0:
        sign = -1
    if sign == -1:
        for i in range(len(all_points)):
            all_points[i] *= -1
        for i in range(len(set1)):
            set1[i] *= -1
        for i in range(len(set2)):
            set2[i] *= -1

    previous_end_point = None

    arrow_color = 'teal'
    for i, point in enumerate(all_points):
        arrow_distance = arrow_distance_orig
        if all_points[i] - offset * sign in set1:
            perfect_alignment = True
        else:
            perfect_alignment = False

        if not perfect_alignment:
            insertions_and_deletions_key = abs(all_points[i] - offset * sign)
            n_insertions = len(insertions_and_deletions[insertions_and_deletions_key][0])
            n_deletions = len(insertions_and_deletions[insertions_and_deletions_key][1])
            arrow_distance += (n_deletions - n_insertions)

        end_point = point+arrow_distance*sign

        tail_length = arrow_distance * .6
        tail_width = .1 
        head_length = (arrow_distance * .4)
        head_width = .2

        if sign == 1:
            bottom_left = point
        else:
            bottom_left = point - tail_length

        tail = Rectangle((bottom_left, y_index - tail_width / 2), tail_length, tail_width, color=arrow_color)
        ax.add_patch(tail)

        if sign == -1 and point < global_min_x:
            global_min_x = point
        elif point > global_max_x:
            global_max_x = point

        triangle = Polygon([[point + (tail_length+head_length) * sign, y_index],
                           [point+(tail_length) * sign, y_index-head_width/2],
                           [point+(tail_length) * sign, y_index+head_width/2]],
                           closed=True, color=arrow_color)
        ax.add_patch(triangle)


        if not perfect_alignment:
            insertions_and_deletions_key = abs(all_points[i] - offset * sign)
            insertions = insertions_and_deletions[insertions_and_deletions_key][0]
            deletions = insertions_and_deletions[insertions_and_deletions_key][1]
            for i in insertions:
                #plt.plot([point+i*sign, point+i*sign], [y_index-.1, y_index+.1], color='blue', lw=1)  # Draw vertical line at midpoint
                ax.plot([point+i*sign, point+i*sign], [y_index-.1, y_index+.1], color='blue', linestyle = '-', lw=1)
            for i in deletions:
                #plt.plot([point+i*sign, point+i*sign], [y_index-.1, y_index+.1], color='gold', lw=1)  # Draw vertical line at midpoint
                ax.plot([point+i*sign, point+i*sign], [y_index-.1, y_index+.1], color='gold', linestyle = '-', lw=1)

        #if previous_end_point is not None and abs(point - previous_end_point) > 1:
        if previous_end_point is not None and abs(point - previous_end_point) > offset_offset+1:
            #ax.plot([previous_end_point, point], [y_index, y_index], linewidth= 3, color='red', linestyle='-')
            ax.plot([previous_end_point+offset_offset/2*sign, point-offset_offset/2*sign], [y_index, y_index], linewidth= 2, color='red', linestyle='-')
        previous_end_point = end_point
        offset += offset_offset


# this is a temporary function that shouldn't be necessary if you filter stuff out as your adding it instead
def final_filter(my_dict, my_filter_options):
    last_diff = sys.maxsize
    last_key = ""
    last_n_subStr = 0
    last_n_matches = 0
    my_arr = list(my_dict.keys())
    required_substring_length = my_filter_options.required_substring_length
    for key in my_arr:
        tmp = my_dict[key]
        if my_filter_options.only_non_overlapping == True and tmp.is_overlap == True:
            del(my_dict[key])
            continue
        if tmp.n_subStr < my_filter_options.min_substring_length:
            del(my_dict[key])
            continue
        if len(tmp.indexes) < my_filter_options.min_matches:
            del(my_dict[key])
            continue
        tmp_diff = tmp.min_gap - tmp.n_subStr
        if required_substring_length > 0:
            if tmp.n_subStr == required_substring_length:
                last_key = key
        elif tmp.min_gap - tmp.n_subStr <= last_diff and tmp.n_subStr > last_n_subStr and len(tmp.indexes) > last_n_matches/2:
            last_diff = tmp_diff
            last_key = key
            last_n_subStr = tmp.n_subStr
            last_n_matches = len(tmp.indexes)
    if required_substring_length > 0:
        if last_key != "":
            my_dict = {last_key: my_dict[last_key]}
            return(my_dict)
        else:
            return None
    elif my_filter_options.is_best_only:
        #my_dict = {last_key: my_dict[last_key]}
        if last_key != "":
            my_dict = {last_key: my_dict[last_key]}
    return(my_dict)

def consolidate(input_s, my_dict, coverage_dict):
    for initial_key in coverage_dict.keys():
        secondary_keys = list(coverage_dict[initial_key].keys())
        secondary_keys.sort()
        for i in range(len(secondary_keys)):
            start = secondary_keys[i]
            end = coverage_dict[initial_key][start]
            if input_s[start:end+1] not in my_dict:
                continue
            for j in range(i+1, len(secondary_keys)):
                cmp_start = secondary_keys[j]
                if cmp_start > end:
                    break
                cmp_end = coverage_dict[initial_key][secondary_keys[j]]
                del my_dict[input_s[cmp_start:cmp_end+1]]
    
    for key in my_dict:
        print(my_dict[key])
        print(key)
        print("\n")
    return

    my_filter_options = filter_options()
    new_dict = final_filter(my_dict, my_filter_options)
    if new_dict == None:
        print("ERROR, no string of the required size was found in this chr end\n")
        return
    output_dict = {}
    for key in new_dict:
        doesnt_fit = True
        #new_dict[key].extra_alignment_indexes, all_extra_alignment_scores = alignment(input_s, key)
        new_dict[key].extra_alignment_indexes, all_extra_alignment_scores, new_dict[key].extra_alignment_insertions_and_deletions  = alignment(input_s, key, my_dict[key])
        new_dict[key].n_extra_alignment = len(new_dict[key].extra_alignment_indexes)
        tmp = new_dict[key]
        if tmp.n_subStr not in output_dict:
            output_dict[tmp.n_subStr] = [[len(tmp.indexes), (len(tmp.indexes)+tmp.n_extra_alignment), [tmp.indexes[0],tmp.indexes[1]]]]
            print(f"{new_dict[key]}all alignment scores: {all_extra_alignment_scores}\n{key}\n")
        else:
            output_dict_val = output_dict[tmp.n_subStr]
            for val in output_dict_val:
                if (val[2][0]<=tmp.indexes[0] and val[2][0]+tmp.n_subStr >= tmp.indexes[0]) or (val[2][1]<=tmp.indexes[0] and val[2][1]+tmp.n_subStr >= tmp.indexes[0]):
                    doesnt_fit = False
                    if val[1] < len(tmp.indexes)+tmp.n_extra_alignment:
                        val[1] = len(tmp.indexes)+tmp.n_extra_alignment
                        print(f"{new_dict[key]}all alignment scores: {all_extra_alignment_scores}\n{key}\n")
                    break
            if(doesnt_fit):
                output_dict_val.append([len(tmp.indexes), (len(tmp.indexes)+tmp.n_extra_alignment), [tmp.indexes[0], tmp.indexes[1]]])
                print(f"{new_dict[key]}all alignment scores: {all_extra_alignment_scores}\n{key}\n")
    if my_filter_options.required_substring_length != -1:
        graph_output(new_dict)

        #alignment(input_s, key)


def expand_dict(input_s, my_dict, coverage_dict):
    input_s_n = len(input_s)
    queue = deque()
    for key in my_dict.keys():
        queue.append(key)
    while queue:
        curr = queue.popleft()
        n_string = my_dict[curr].n_subStr
        indexes = my_dict[curr].indexes
        n_indexes = len(indexes)
        s = ""
        tmp_dict = {}
        for index in indexes:
            r_index = index + n_string
            if r_index >= input_s_n:
                continue
            s = input_s[index:r_index+1]
            if s in my_dict:
                gap = index - my_dict[s].indexes[-1]
                if gap < my_dict[s].min_gap:
                    my_dict[s].min_gap = gap
                    my_dict[s].is_overlap = gap < len(s)
                my_dict[s].indexes.append(index)
                tmp_dict[s] = index
            elif s in tmp_dict:
                first_val = tmp_dict[s]
                gap = index - first_val
                is_overlap = gap < len(s)
                my_dict[s] = self_search_type(len(s), [first_val, index], gap, is_overlap)
                queue.append(s)
            else:
                tmp_dict[s] = index
        for key in tmp_dict.keys():
            if key in my_dict:
                tmp = my_dict[key]
                n_tmp_indexes = len(tmp.indexes)
                start_index = tmp.indexes[0]
                if n_tmp_indexes not in coverage_dict:
                    coverage_dict[n_tmp_indexes] = {start_index: start_index + tmp.n_subStr - 1}
                else:
                    if  start_index in coverage_dict[n_tmp_indexes]:
                        end_index = coverage_dict[n_tmp_indexes][start_index]
                        sub_str = input_s[start_index:end_index+1]
                        if my_dict[sub_str].is_overlap == my_dict[key].is_overlap:
                            del my_dict[input_s[start_index:end_index+1]]
                    coverage_dict[n_tmp_indexes][start_index] = start_index + tmp.n_subStr - 1
    consolidate(input_s, my_dict, coverage_dict)

            

def self_search(input_s, min_length = 50):
    #global global_counter
    #global_counter += 1
    #if global_counter%2 == 0:
        #x = "L"
    #else:
        #x = "R"
    #print(f"<<<<<<<<<<<CHR{global_counter//2+1}{x}>>>>>>>>>>>>>>>>")
    #if input_s == "":
        #print("skipping, emptry string was given\n")
        #return
    my_dict = {}
    coverage_dict = defaultdict(dict)
    n = len(input_s)
    tmp_dict = defaultdict(list)
    for i in range(0, n-min_length+1):
        sub_str = input_s[i:i+min_length]
        tmp_dict[sub_str].append(i)
    # coverage_dict is a dictionary of dictionaries. First value refers to how many matches for a particular value, second value is the start of the coverage, and the dictionary value is the end of the coverage
    for key in tmp_dict.keys():
        tmp_arr = tmp_dict[key]
        len_arr = len(tmp_arr)
        if len_arr > 2:
            min_gap = sys.maxsize
            for i in range(1, len_arr):
                min_gap = min(min_gap, tmp_arr[i] - tmp_arr[i-1])
            first_index = tmp_arr[0]
            coverage_dict[len_arr][first_index] = first_index+min_length - 1
            overlap = min_gap < min_length
            my_dict[key] = self_search_type(min_length, tmp_arr, min_gap, overlap)
    if len(list(my_dict.keys())) == 0:
        print("no loops of the minimum length were found\n")
        return
    #for key in my_dict:
        #print(my_dict[key])
        #print(key)
        #print("\n")
    #return
    expand_dict(input_s, my_dict, coverage_dict)


def read_fasta(file_path):
    fasta_dict = {}
    with open(file_path, 'r') as f:
        identifier = ''
        sequence_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if identifier:
                    fasta_dict[identifier] = ''.join(sequence_lines)
                    sequence_lines = []
                identifier = line[1:]
            else:
                sequence_lines.append(line)
        if identifier:
            fasta_dict[identifier] = ''.join(sequence_lines)
    return fasta_dict

def pre_processing(input_s, dict_key=-1, file_path = "6991_only_telomeres_new.fasta"):
    global global_counter
    global_counter += 1
    if global_counter%2 == 0:
        x = "L"
    else:
        x = "R"
    print(f"<<<<<<<<<<<CHR{global_counter//2+1}{x}>>>>>>>>>>>")
    if input_s == "":
        print("skipping, emptry string was given\n")
        return
    dict = read_fasta(file_path) 
    if dict_key == -1:
        num = (global_counter) // 2 + 1
        dict_key = f"chr{num}{x}_Telomere_Repeat"
    known_value = dict[dict_key]

    if input_s[0] == "A" or input_s[0] == "C":
        print("converting string to its reverse complement")
        tmp = Seq(input_s)
        input_s = str(tmp.reverse_complement())

    input_n = len(input_s)
    known_value_n = len(known_value)

    output = 0
    tmp_value = 0

    for i in range(5):
        for j in range(min(known_value_n, input_n)):
            if known_value[j] == input_s[j+i]:
                tmp_value += 1
            else:
                break
        output = max(output, tmp_value)
        if output > 5:
            #return(output+i)
            print(f"deleting first {output+i} characters")
            self_search(input_s[output+i:])
            return
        tmp_value = 0
    self_search(input_s)


def main():


    self_search("ACCCACACCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACCACACCCACACCCACACACCACACCCACACCCACACACACCACACCCACACCCACACCCACACACCACACCCACACCCACACACCACACCCACACCCACACCCACACACCCACACCACACCCACACACCACACCCACACCACACCCACACACCCACACACACACACCCACACACCACACCCACACCACACCCACACCACACCCACACCCACACCCACACCACACCCACACACACCACACCCACACCCACACACCACACACAC")
    return

    graph_setup()

    #IT156

    #chr1L
    pre_processing("ACCCACACCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACCACACCCACACCCACACACCACACCCACACCCACACACACCACACCCACACCCACACCCACACACCACACCCACACCCACACACCACACCCACACCCACACCCACACACCCACACCACACCCACACACCACACCCACACCACACCCACACACCCACACACACACACCCACACACCACACCCACACCACACCCACACCACACCCACACCCACACCCACACCACACCCACACACACCACACCCACACCCACACACCACACACAC")

    #chr1R
    pre_processing("GGGTGTGTGTGTGTGGGTGTGGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTAGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGGTGTGTGGGTGTGGTGTGGTGTGGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGGGTGGGGGTGGTGTGGGGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGTGTGGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGT")

    #chr2L
    pre_processing("TGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGGGTGTGTGTGTGTGGGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGGGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGGGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGGGTGTGCGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGGGGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGT")

    #chr2R
    pre_processing("ACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACCCACACACCACACCCACACCCACACACCACACCCACAC")

    #chr3L
    pre_processing("GGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTG")

    #chr3R
    pre_processing("GTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGCGGTTGTGTGTGTGGGTGTGGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTG")

    #chr4L
    pre_processing("GGGTGTGGGTGTGTGTGTGTGGTGTGTGTGTGGTGTGGGTGTGGGTGTGGGGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGTGTGGGGGTGTGTGTGGGTGTGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGGGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGTGGGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGGTGTGGGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGT")

    #chr4R
    pre_processing("GTGTGGTGTGTGGGTGTGGGTGTGTGGTGTGTGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGGGTGTGGTGGGTGTGGTGTGGTGTGTGTGGGTGTGGTGTGGTGTGT")

    #chr5L
    pre_processing("GTGTGTGGTGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGTGGTGTGGTGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGGGTG")

    #chr5R
    pre_processing("GTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGT")

    #chr6L
    pre_processing("GTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGTGTGTGTGTGGGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTG")

    #chr6R
    pre_processing("")
    #chr7L
    pre_processing("")
    #chr7R
    pre_processing("")
    #chr8L
    pre_processing("GTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGT")
    #chr8R
    pre_processing("GTGTGGGTGTGGTGTGTGGGTGTGGGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTG")
    #chr9L
    pre_processing("CACACACCCACACCACACCACACCACACCCACACACACACCCACACACACACCACACCCACACCCACACCCACACACACACACACCACACCCACACCACACCACACCCACACACACACCCACACCACACACCCACACCACACCACACCCACACACACACCCACACACACACACCCACACCACACCCACACACACACCACACACACACACCCACACCACACCCCACCCACACACCCACACCACACCCACACCACACCACACCCACACACACACACACCCACACACACCCACACCCCACCCACACCACACCACACCCACACACACACGCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCCACACACACACACCCACACCACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCCACACCACACCCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACACACCACACCCACACACCCACACACCACACACCACACAC")
    #chr9R
    pre_processing("GTGTGTGGGTGTGGGTGGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGTGGTGTGGGTGTGGTGTGGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGT")
    #chr10L
    pre_processing("")
    #chr10R
    pre_processing("GTGTGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGGGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGG")
    #chr11L
    pre_processing("CACCACACCCACACACACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACCCACACACCACACCCACACCCACACACCACACCCAC")
    #chr11R
    pre_processing("CACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACACCCACACACACACACACCACACCCACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACACACCACACCCACACACACACACCCACACACACACACCCACACACACACACCCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACACCCACACACACACACACCACACCCACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACCACACCCACACACCCACAACCC")
    #chr12L
    pre_processing("ACCCACACCACACCCACCACACCACACCCACACACACACCCACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCCACACACACACACCCACACCACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCCACACCACACCCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACACACCACACCCACACACCCACACACCACACACCACACAC")
    #chr12R
    pre_processing("")

    #chr13L
    pre_processing("")
    #chr13R
    pre_processing("GGGTGTGGGTGTGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTG")
    #chr14L
    pre_processing("TGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGTGTGGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGTGGGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGGTGTGTGTGTGTGTGGTATGTGTGTGTGTGGGTGTGGTGGTGTGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGTGGTGTGGTGTGGTGGTGTGGTGGGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGAGTGTGTGTGGTGTGGTGTGTGTGGTGTGGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGGTGTGGTGGGGTGGTGGGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGTGTGGGTGTGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGGGTGTGGTGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGGGTGTGTGTGTGGGTGTGGTG")
    #chr14R
    pre_processing("")
    #chr15L
    pre_processing("")
    #chr15R
    pre_processing("GTGTGTGGGTGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGGTGTGGGTGTGTGTGGGTGTGGGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGGTGTGGGTGTGTGTGTGGTGTGTGTGTGTGGTGTGGTGTGTGTGTGGGTGTGGTGTGGGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGTGTGGTGTGTGTGGGTGTGGTGTGGTGTGTGTGTGGTGTGTGTGTGTGGGTGTGTGTGTGGTGTGGGTGTGGGTGTGGGTGTGTGGTGTGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGGGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGGGTGTGTGGTGTGGGTGTGGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGGGGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGGGGTGTGGGTGTGTGTGTGTGTGGGTGTGGTGGGTGTGGTGTGGGTGTGGGGGTGGTGTGTGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGTGTG")
    #chr16L
    pre_processing("TGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGT")
    #chr16R
    pre_processing("GTGTGGGTGTGGTGGTGTGGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGGGTGTGT")

    # Set the figure’s DPI to match the interactive display
    #plt.gcf().set_dpi(plt.gcf().get_dpi())  # This will use the current DPI

    # Save the plot with tight bounding box to capture layout
    #plt.savefig("my_plot.png", dpi=plt.gcf().get_dpi(), bbox_inches='tight')
    ax.set_xlim(global_min_x, global_max_x)

    fig.savefig(global_figure_name, dpi=300, bbox_inches='tight')

    #plt.show()

if __name__ == "__main__":
    main()