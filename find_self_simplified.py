from collections import deque, defaultdict
import sys
from Bio import Align
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, Polygon

fig, ax = plt.subplots(figsize=(16,4))

global_counter = -1 

class filter_options:
    def __init__(self, only_non_overlapping = True, min_substring_length = 50, min_matches = 3, is_best_only = True):
        self.only_non_overlapping = only_non_overlapping
        self.min_substring_length = min_substring_length
        self.min_matches = min_matches
        self.is_best_only = is_best_only

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
        avg_dist = self.average_distance()
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
            score = alignments[0].score
            # below if statement determines what qualifies as a good alignment score
            #if score >= key_n - key_n/20:
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

                new_entry_start = str_dict[str][0] + alignments[0].coordinates[0][0] + key_n + 1
                new_entry_end = str_dict[str][1]
                new_str = input_s[new_entry_start:new_entry_end]
                str_dict[new_str] = [new_entry_start, new_entry_end]
                queue.append(new_str)

                coordinates = alignments[0].coordinates
                insertions = 0
                insertions_indexes = []
                deletions=0
                deletions_indexes = []

                # Iterate over each alignment block in coordinates
                for i in range(1, len(coordinates[0])-2, 2):
                    # Extract the end of the current block and the start of the next
                    seq1_end, seq1_next = coordinates[0][i], coordinates[0][i+1]
                    seq2_end, seq2_next = coordinates[1][i], coordinates[1][i+1]

                    # Count gaps between alignment blocks
                    if seq1_end != seq1_next:  # Deletion in seq2
                        deletions += seq1_next - seq1_end
                        for i in range(seq1_next-seq1_end):
                           deletions_indexes.append(seq2_end + i)
                    if seq2_end != seq2_next:  # Insertion in seq1
                        insertions += seq2_next - seq2_end
                        for i in range(seq2_next-seq2_end):
                           insertions_indexes.append(seq2_end + i)

                dict_insertions_and_deletions[str_dict[str][0]+alignment_start] = [insertions_indexes, deletions_indexes]

            arr.append(score)
    return(good_alignment_starting_pos, arr, dict_insertions_and_deletions)

    #print(f"number of good extra alignments: {counter}\nextra alignment scores: {arr}\ngood alignment starts: {good_alignment_starting_pos}\n")
    #for alignment in good_alignment_arr:
        #print(alignment)
    #print(f"number of good extra alignments: {counter}, extra alignment scores: {arr}\nThe best extra Alignment that was found:\n{best_alignment}\n")

def graph_setup():
    plt.figure(figsize=(8,6))
    # Set a fixed DPI for both display and saving
    fig = plt.gcf()
    #plt.figure(dpi=150)
    #display_dpi = fig.get_dpi()  # Get the current display DPI
    #fig.set_dpi(display_dpi)
    #plt.title(0, 8.5, "IT148")
    plt.text(0, 8, "IT148", ha='center', va='center', fontsize=12, color='black')
    plt.text(0, 7, "chr1", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 6.5, "chr2", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 6, "chr3", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 5.5, "chr4", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 5, "chr5", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 4.5, "chr6", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 4, "chr7", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 3.5, "chr8", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 3, "chr9", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 2.5, "chr10", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 2, "chr11", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 1.5, "chr12", ha='center', va='center', fontsize=10, color='black')
    plt.text(0, 1, "chr13", ha='center', va='center', fontsize=10, color='black')

    ax.text(0, 8, "IT148", ha='center', va='center', fontsize=12, color='black')
    ax.text(0, 7, "chr1", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 6.5, "chr2", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 6, "chr3", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 5.5, "chr4", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 5, "chr5", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 4.5, "chr6", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 4, "chr7", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 3.5, "chr8", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 3, "chr9", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 2.5, "chr10", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 2, "chr11", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 1.5, "chr12", ha='center', va='center', fontsize=10, color='black')
    ax.text(0, 1, "chr13", ha='center', va='center', fontsize=10, color='black')

    # Add x-axis label and set limits
    plt.xlabel("X-axis")
    ax.set_xlabel("X-axis")
    plt.ylim(.5, 7.5)  # Set y-axis limits to keep the line flat
    ax.set_ylim(.5, 7.5)
    ax.set_xlim(-6000, 5000)


    # Remove y-axis and add legend
    plt.gca().get_yaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

def graph_output(my_dict):
    offset = 300 
    global global_counter
    y_index = (14 - (global_counter//2)) / 2
    key = list(my_dict.keys())[0]
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
            all_points.append(set1_tmp + offset)
            set1_ptr += 1
        else:
            all_points.append(set2_tmp + offset)
            set2_ptr += 1
    if set1_ptr == set1_n:
        for i in range(set2_ptr, set2_n):
            all_points.append(set2[i] + offset)
    else:
        for i in range(set1_ptr, set1_n, 1):
            all_points.append(set1[i] + offset)
    
    arrow_style = "->"
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

    arrow_distance_orig = my_dict[key].n_subStr
    #arrow_distance = my_dict[key].n_subStr

    previous_end_point = None

    # Create the 1D plot along the x-axis
    #plt.plot(all_points, [y_index] * len(all_points), 'o', color='black', markersize=0, label='Data Points')  # Black dots for all data points

    # Add arrows between each point, with color based on originating set
    arrow_color = 'teal'
    for i, point in enumerate(all_points):
        arrow_distance = arrow_distance_orig
        if all_points[i] - offset * sign in set1:
            #arrow_color = 'teal'
            #perfect_alignment = 'teal'
            perfect_alignment = True
        else:
            #arrow_color = 'skyblue' 
            #perfect_alignment = 'skyblue' 
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

        #arrowhead = FancyArrowPatch((point + tail_length, y_index), (end_point, y_index),
                            #color="blue", arrowstyle="->", mutation_scale=100, lw=2)
        #ax.add_patch(arrowhead)

        triangle = Polygon([[point + (tail_length+head_length)*.9 * sign, y_index],
                           #[point+(tail_length+head_length) * sign, y_index-head_width/2],
                           [point+(tail_length) * sign, y_index-head_width/2],
                           #[point+(tail_length+head_length) * sign, y_index+head_width/2]],
                           [point+(tail_length) * sign, y_index+head_width/2]],
                           closed=True, color=arrow_color)
        ax.add_patch(triangle)

        #arrow = FancyArrowPatch((end_point, y_index), (point, y_index),
                        #color=arrow_color, linewidth=3, arrowstyle="->",
                        #mutation_scale=10)  # Adjust mutation_scale for arrowhead

        #ax.add_patch(arrow)

        plt.annotate('', xy=(end_point, y_index), xytext=(point, y_index),
                    arrowprops=dict(arrowstyle=arrow_style, color=arrow_color, lw=1.5,
                                    mutation_scale=10, shrinkA=0, shrinkB=0))


        if not perfect_alignment:
            insertions_and_deletions_key = abs(all_points[i] - offset * sign)
            insertions = insertions_and_deletions[insertions_and_deletions_key][0]
            deletions = insertions_and_deletions[insertions_and_deletions_key][1]
            for i in insertions:
                plt.plot([point+i*sign, point+i*sign], [y_index-.1, y_index+.1], color='blue', lw=1)  # Draw vertical line at midpoint
                ax.plot([point+i*sign, point+i*sign], [y_index-.1, y_index+.1], color='blue', linestyle = '-')
            for i in deletions:
                plt.plot([point+i*sign, point+i*sign], [y_index-.1, y_index+.1], color='gold', lw=1)  # Draw vertical line at midpoint
                ax.plot([point+i*sign, point+i*sign], [y_index-.1, y_index+.1], color='gold', linestyle = '-')

        if previous_end_point is not None:
            plt.plot([previous_end_point, point], [y_index, y_index], color='red', lw=1)
            ax.plot([previous_end_point, point], [y_index, y_index], linewidth= 3, color='red', linestyle='-')
        previous_end_point = end_point


# this is a temporary function that shouldn't be necessary if you filter stuff out as your adding it instead
def final_filter(my_dict, my_filter_options):
    last_diff = sys.maxsize
    last_key = ""
    last_n_subStr = 0
    last_n_matches = 0
    my_arr = list(my_dict.keys())
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
        if tmp.min_gap - tmp.n_subStr <= last_diff and tmp.n_subStr > last_n_subStr and len(tmp.indexes) > last_n_matches/2:
            last_diff = tmp_diff
            last_key = key
            last_n_subStr = tmp.n_subStr
            last_n_matches = len(tmp.indexes)
    if my_filter_options.is_best_only:
        #my_dict = {last_key: my_dict[last_key]}
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

    my_filter_options = filter_options()
    new_dict = final_filter(my_dict, my_filter_options)
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
    global global_counter
    global_counter += 1
    if input_s == "":
        return
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
        return
    if global_counter%2 == 0:
        x = "L"
    else:
        x = "R"
    print(f"<<<<<<<<<<<CHR{global_counter//2+1}{x}>>>>>>>>>>>>>>>>")
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

def pre_processing(input_s, dict_key, file_path = "6991_only_telomeres_new.fasta"):
    dict = read_fasta(file_path) 
    known_value = dict[dict_key]

    if input_s[0] == "A" or input_s[0] == "C":
        print("reversing string")
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
    graph_setup()

    #48 - chr1l - has cycles
    #self_search("GTGTGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGGGTGGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGGTGTGGTGTGTGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGTGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGGGTGTGGGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGTGTGGGGTGTGGTGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGGGGTGTGGTGTGTGTGTGTGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTG", "chr1L_Telomere_Repeat")
    self_search("GTGTGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGGGTGGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGGTGTGGTGTGTGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGTGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGGGTGTGGGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGTGTGGGGTGTGGTGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGGGGTGTGGTGTGTGTGTGTGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTG")

    #48 - chr1r - has cycles
    self_search("GGTGTGTGTGGGTGTGGTGTGGGTGTGGTTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGTGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGGTGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGGTGTGTGTGTGTGGGTGTGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTG")

    #48 - chr2l - has cycles
    self_search("ACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACCCACACACCACACCCACACACCCACACACCACACCCACACACCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACCCCACACCCACACACCACACCCACACACCCACACACCCACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCACACCCACACACACACCACACCCCACACACACACCACACCCACACACCCACACCCACACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACACCACACACACCACACCCACACACACACACACCCACACACCCACACACCACACCCACACACACACACACACCACACACACACACACACACACACCACACACCCACACACCACACCCACACACACCCACACCACACACACACACACCACACACACACACACCCACACACCACACCCACACACCACACCCACACACCCACACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACACCACACACACCACACCCACACACCCACACCCACACCCACACACACACCCCACACACCACACCCACACACCCACACCCACACCACACCCACACACACCACACACCACACCCACACACCACACCCACACACCACCCACACCCACACAC")
    #pre_processing("ACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACCCACACACCACACCCACACACCCACACACCACACCCACACACCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACCCCACACCCACACACCACACCCACACACCCACACACCCACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCACACCCACACACACACCACACCCCACACACACACCACACCCACACACCCACACCCACACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACACCACACACACCACACCCACACACACACACACCCACACACCCACACACCACACCCACACACACACACACACCACACACACACACACACACACACCACACACCCACACACCACACCCACACACACCCACACCACACACACACACACCACACACACACACACCCACACACCACACCCACACACCACACCCACACACCCACACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACACCACACACACCACACCCACACACCCACACCCACACCCACACACACACCCCACACACCACACCCACACACCCACACCCACACCACACCCACACACACCACACACCACACCCACACACCACACCCACACACCACCCACACCCACACAC", "chr2L_Telomere_Repeat")

    #48 - chr2r - has cycles
    self_search("ACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACCCACACACCACACCCACACACCCACACACCACACCCACACACCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACCCCACACCCACACACCACACCCACACACCCACACACCCACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCACACCCACACACACACCACACCCCACACACACACCACACCCACACACCCACACCCACACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACACCACACACACCACACCCACACACACACACACCCACACACCCACACACCACACCCACACACACACACACACCACACACACACACACACACACACCACACACCCACACACCACACCCACACACACCCACACCACACACACACACACCACACACACACACACCCACACACCACACCCACACACCACACCCACACACCCACACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACACCACACACACCACACCCACACACCCACACCCACACCCACACACACACCCCACACACCACACCCACACACCCACACCCACACCACACCCACACACACCACACACCACACCCACACACCACACCCACACACCACCCACACCCACACAC")

    #48 - chr3l (ending looks kinda weird in this one)
    self_search("GGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGT")

    # chr3r - some long cycles with lots of overlap
    self_search("ACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACCCACACACACACCCACACCACACCCACACACCACACCCACACCCACACACCCACACACCACACCACTAACCCTAACACTACCC")

    # chr4l (short telomer) - none
    self_search("")

    #chr4r (no telomer)
    self_search("")

    # chr5L (no telomre)
    self_search("")

    # chr5R (multiple contigs)
    self_search("")

    # chr6L - ??
    self_search("GTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGGGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGGTGGTGTGTGGTGTGGTGTGGGTGTGTGTGTGGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGGTGTGTGTGGTGTGGTGTGTGGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGTGGGGTGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGGTGTGTGGTGGGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTG")
    #self_search("")

    # chr6R - multiple congigs
    self_search("")

    # chr7L - has cycles
    self_search("GTGTGTGGGTGTGGGTGTGGTGTGTGGTGTGGGTGTGTGGGTAGTGTGTGGGTGTGGGTGTGGTGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGTGTG")

    # chr7R - some
    self_search("GTGTGTGGGTGTGGGTGTGGGTGTGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGGGTGTGGT")
    #self_search("GTGTGTGGGTGTGGGTGTGGGTGTGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGGGTGTGGT")

    # chr8L - has cycles
    self_search("ACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCACACACACACCACACCCACACCACACCCACACACCCACACACACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCCACCACACACCACACCCACCCACACACCACACCCACACACCCACACACCACACCCCACACCCACACCCACACACACACCACACCCCACACACCACACCCCACCCACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACACACACACACACACCACACCCCACACCCCACACACCACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCACACACACCACACCCACACCCACACACCCACACCCCACACCCACACCACACCCACACCCACACACCCACACCCACACACCCACACACCACACCCCACACACCACACCCACACACCACACCCACACCCACACACACACCACACCCACACACCACACACCACACCCACACCACACCCACACACCCACACACCACACCACACACACCACACCCACACACCCACACACCACACCCACACACACACCACACCCACACACCACACCCACACCCACACCACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCCACACCCACACACCACACCACACCACACCCACACACCACACCCACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACCACACACCCACACACCACACCACACCACACACCACACACCACACCCACACCCACACACCACACCCACACACCCACACCCACACCCACACAC")

    # chr8R - some
    #self_search("GTGTGGGTGTGGTGGTGTGGGTGTGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTG")
    self_search("GTGTGGGTGTGGTGGTGTGGGTGTGTGTGGTGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTG")

    # chr9L emptry
    self_search("")

    # chr9R - has cycles
    self_search("GTTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGTGTGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTG")

    # chr10L - has cycles
    self_search("ACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACACACCCACACACCCACACCCACACCCACACACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCCACACACCACACCCACACACCCACACACCACACCACACACACACCACACCCACACACACCACACCCACACACCACACACACCACACCCACACACCCACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACCCACACCCACACACCACACCAC")

    # chr10R
    self_search("GTTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGTGTGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGGGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGTGTGTGGGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGG")

    #chr11L - has cycles
    self_search("GGGTTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTG")

    #chr11R - none
    self_search("")

    #chr12L - has cycles
    self_search("CCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACACCACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCCACCACCACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACCACACCCACACACCCACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACCCACACCCACACCACACCCACACACCACACCCACACCCCACACACCACACCCACACACCCACACCCACACCCACACCACACACCACACCCCACACACACACCACACCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACCACACACACACCACACCCCACACACACACCACACCCACCCACACCCACACCCCACACACCACACCCACACACCACACCCACACCCCACACACCACACCCCACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACCACACCCACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCCACACCCACCACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACACCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACACACACCACACACACACCACACCCACACCACACCCACACCCACACCCACACACCCACACCCACACAC")

    #Chr12R super short telomer, 2 contigs
    self_search("")

    #chr13L - has cycles
    self_search("CCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACCACACCCACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACCCACACACACACCCACACACCCACACCCACACACCCACACCCACACACCACACACCCACACACCACACCACACCCACACCACACCCACCACACCCACACACCACACCCACACCACACCCACACACCCACACCCACACCCACACACCACACCCACACACCCACACACCACACCACACCACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACACACACACCACACCCACACACCCACACCCACACCCACACCCACACACCACACCCACACACCCACACACCACACCCACAC")

    #chr13R - has cycles
    self_search("CACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACACACACCACACCCACACACACCACACCCACACACACACCACACCCACACACACACACCACACCCACACCCTAACAC")

    ## IT145:

    #chr1L - weird double telomer, cycles found
    #self_search("GTGTTGAGGTAGTGTTAGTGTGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGGGTGTGTGTGGTGTGGTGTGGTGTGGGTGTGGTGTGGGTGTGGGTGTGGGTGTGGTGTGGGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGTGGTGTGGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGTGGATGTGGGTGTGGTGTGGGTGTGGTGTGGTGTGGAGCAATACGTATACTGTGACGATGTACTTCGCTACAGTATACGCATTGCTCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACCACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACACCACACCCACATCCACACACCACACCCACACACCCACACACCACACCCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCCACACACCACACCCACACCACACCCCACACCACACATCCACACCACCCACACACCCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCAC")

    #chr2L
    #self_search("ACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCCACACCACACCCACATCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACATCCACACACCCACACACCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACACCCACACACCACACCCACACCACACCCACACCACACCCACATCCACACACCCACACACCCACACACCCACACCCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCACACCCACACCCACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCACCCACACCCACACAC")

    # Set the figure’s DPI to match the interactive display
    #plt.gcf().set_dpi(plt.gcf().get_dpi())  # This will use the current DPI

    # Save the plot with tight bounding box to capture layout
    #plt.savefig("my_plot.png", dpi=plt.gcf().get_dpi(), bbox_inches='tight')
    fig.savefig('my_plot.png', dpi=100, bbox_inches='tight')

    plt.show()

if __name__ == "__main__":
    main()