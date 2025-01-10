from collections import deque, defaultdict
import sys
import copy
import time
import re
import argparse
import os
from Bio import Align, SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
from contextlib import redirect_stdout


def is_circular_rearrangement(s1, s2):
    if len(s1) != len(s2):
        return False

    s2_double = s2 + s2

    if s1 in s2_double:
        return True
    return s1[::-1] in s2_double


def convert_num_to_chr_end(num):
    tmp = (num+2)//2
    if num % 2 == 0:
        return f"{tmp}L"
    else:
        return f"{tmp}R"


class self_search_type:
    def __init__(self, n_subStr, indexes, min_gap, has_overlap, offset):
        self.n_subStr = n_subStr
        self.indexes = indexes
        self.min_gap = min_gap
        self.has_overlap = has_overlap
        self.offset = offset
        # self.n_extra_alignment = 0
        self.extra_alignment_indexes = []
        self.extra_alignment_insertions_and_deletions = {}
        # self.parent_key = parent_key

    def add_offset(self):
        for i in range(len(self.indexes)):
            self.indexes[i] += self.offset
        for i in range(len(self.extra_alignment_indexes)):
            self.extra_alignment_indexes[i] += self.offset
        # for key in self.extra_alignment_insertions_and_deletions.keys():
            # for i in range(len(self.extra_alignment_insertions_and_deletions[key])):
            # for j in range(len(self.extra_alignment_insertions_and_deletions[key][i])):
            # self.extra_alignment_insertions_and_deletions[key][i][j] += self.offset
            # self.extra_alignment_insertions_and_deletions[key][i] += self.offset
            # for i in range(len(self.extra_alignment_insertions_and_deletions[key][1])):
            # self.extra_alignment_insertions_and_deletions[key][i] += self.offset

    def __repr__(self):
        n = self.n_subStr
        # avg_dist = self.average_distance()
        # return f"substring length={self.n_subStr}, number of matches={len(self.indexes)}, min_gap={self.min_gap}, overlap={self.has_overlap}, number of extra alignments={self.n_extra_alignment}\nstarting indexes: {self.indexes}\nextra alignemtn indexes: {self.extra_alignment_indexes}\n"
        return f"substring length={self.n_subStr}, number of matches={len(self.indexes)}, min_gap={self.min_gap}, overlap={self.has_overlap}, number of extra alignments={len(self.extra_alignment_indexes)}\nstarting indexes: {self.indexes}\nextra alignemtn indexes: {self.extra_alignment_indexes}\n"


class user_input:
    def __init__(self, **kwargs):
        # Dynamically set attributes based on the provided keyword arguments
        self.__dict__.update(kwargs)


class find_loops:

    def __init__(self, input_s, user_input_obj, offset):
        self.input_s = input_s
        self.my_dict = {}
        self.coverage_dict = defaultdict(dict)
        # self.min_length = min_length
        self.user_input_obj = user_input_obj
        self.offset = offset

        self.aligner = Align.PairwiseAligner()

        # Local (Smith-Waterman) alignment; use 'global' for global (Needleman-Wunsch)
        self.aligner.mode = 'local'
        self.aligner.match_score = 1  # Score for a match
        self.aligner.mismatch_score = -(1)  # Penalty for a mismatch
        self.aligner.open_gap_score = -(1)  # Penalty for opening a gap
        self.aligner.extend_gap_score = -(.5)  # Penalty for extending a gap
        # aligner.query_left_open_gap_score = 0
        # aligner.query_right_open_gap_score = 0
        # aligner.target_left_open_gap_score = 0
        # aligner.target_right_open_gap_score = 0

    def expand_dict(self):
        input_s = self.input_s
        my_dict = self.my_dict
        input_s_n = len(input_s)
        queue = deque()
        for key in my_dict.keys():
            queue.append(key)
        while queue:
            curr = queue.popleft()
            # if statement used for immediately removing sequences that contain overlap
            if my_dict[curr].has_overlap:
                del my_dict[curr]
                continue
            n_string = my_dict[curr].n_subStr
            indexes = my_dict[curr].indexes
            # n_indexes = len(indexes)
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
                        my_dict[s].has_overlap = gap < len(s)
                    my_dict[s].indexes.append(index)
                    tmp_dict[s] = index
                elif s in tmp_dict:
                    first_val = tmp_dict[s]
                    gap = index - first_val
                    has_overlap = gap < len(s)
                    my_dict[s] = self_search_type(
                        len(s), [first_val, index], gap, has_overlap, self.offset)
                    queue.append(s)
                else:
                    tmp_dict[s] = index
            for key in tmp_dict.keys():
                if key in my_dict:
                    # if statement used for imediately removing sequences with overlap (added by previous loop)
                    if not my_dict[key].has_overlap:
                        curr_dict_value = my_dict[key]
                        n_curr_dict_value_indexes = len(
                            curr_dict_value.indexes)  # remove this variable
                        if n_curr_dict_value_indexes == len(my_dict[curr].indexes):
                            del my_dict[curr]

    def self_search(self, input_s=""):
        min_length = self.user_input_obj.min_length
        input_s = self.input_s
        n = len(input_s)
        tmp_dict = defaultdict(list)

        # this loops creates tmp_dict which contains all possible substrings of size min_length as the dictionary key, and the indexes of that substring as the value
        for i in range(0, n-min_length+1):
            sub_str = input_s[i:i+min_length]
            tmp_dict[sub_str].append(i)

        # covrage_dict is a dictionary of dictionaries. First value refers to how many matches for a particular value, second value is the start of the coverage, and the dictionary value is the end of the coverag\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\e

        # following code is for converting
        last_index = None
        last_key = None
        should_skip = False
        sequential_englargments_counter = 1

        # for i in range(0, n-min_length+1):
        for key in tmp_dict:
            # sub_str = input_s[i:i+min_length]
            sub_str = key
            tmp_arr = tmp_dict[sub_str]
            len_arr = len(tmp_arr)

            # if the array from our tmp_dict value only contains 1 value then that means it didn't repeat at all and we are not interested in it.
            # following if conditional also determines the minimum number of matches that must be found in order to continue
            if len_arr > 2:

                # if last_index != None and last_index == counter-1:
                if last_key:
                    prev_arr = self.my_dict[last_key].indexes
                    len_prev = len(prev_arr)
                    has_all_same_locations = True
                    should_skip = False

                    # if len_arr <= len_prev:
                    if len_arr == len_prev:
                        curr_arr_ptr = 0
                        prev_arr_ptr = 0
                        while curr_arr_ptr < len_arr:
                            if prev_arr_ptr >= len_prev:
                                has_all_same_locations = False
                                break

                            else:
                                sum = prev_arr[prev_arr_ptr] + \
                                    sequential_englargments_counter
                                if sum > tmp_arr[curr_arr_ptr]:
                                    has_all_same_locations = False
                                    break

                                elif sum == tmp_arr[curr_arr_ptr]:
                                    curr_arr_ptr += 1

                            prev_arr_ptr += 1

                        if has_all_same_locations:
                            if len_arr != len_prev:
                                sequential_englargments_counter += 1
                                continue
                            should_skip = True
                            # tmp_object = self.my_dict[last_key]
                            tmp_object = copy.deepcopy(self.my_dict[last_key])
                            had_overlapping = tmp_object.min_gap < tmp_object.n_subStr
                            tmp_object.n_subStr += 1
                            tmp_object.has_overlap = tmp_object.min_gap < tmp_object.n_subStr
                            new_key = str(last_key)+sub_str[-1]
                            self.my_dict[new_key] = tmp_object
                            if had_overlapping == tmp_object.has_overlap:
                                del self.my_dict[last_key]
                            last_key = new_key
                            sequential_englargments_counter += 1

                # following if statement is for converting the basic array that serves as the value to our tmp_dict dictionary, to a more robust object type that can tonain more information, and we store this object type in the my_dict dictionary
                if not should_skip:
                    min_gap = sys.maxsize
                    for j in range(1, len_arr):
                        min_gap = min(min_gap, tmp_arr[j] - tmp_arr[j-1])
                    first_index = tmp_arr[0]
                    self.coverage_dict[len_arr][first_index] = first_index + \
                        min_length - 1
                    overlap = min_gap < min_length
                    self.my_dict[sub_str] = self_search_type(
                        min_length, tmp_arr, min_gap, overlap, self.offset)
                    last_key = sub_str
                    sequential_englargments_counter = 1

        if len(list(self.my_dict.keys())) == 0:
            # print("no loops of the minimum length were found\n")
            return

    def alignment(self, key):
        my_dict_entry = self.my_dict[key]
        alignments = my_dict_entry.indexes
        key_n = len(key)

        last_val = 0
        str_dict = {}
        queue = deque()
        for alignment in alignments:
            tmp = self.input_s[last_val:alignment]
            str_dict[tmp] = [last_val, alignment]
            queue.append(tmp)
            last_val = alignment+key_n

        # if last_val < len(self.input_s):
            # tmp = self.input_s[last_val:]
            # str_dict[tmp] = [last_val, ""]
            # queue.append(tmp)

        good_alignment_arr = []
        good_alignment_starting_pos = []
        dict_insertions_and_deletions = {}
        arr = []
        counter = 0
        while queue:
            str = queue.popleft()
            if str != "":
                alignments = self.aligner.align(str, key)
                # try except block is necessary cause of overflow if there are too many optimal alignments
                try:
                    n_alignments = len(alignments)
                except OverflowError:
                    n_alignments = 1

                # if statement is necessary because sometimes there arent any good alignments
                if n_alignments > 0:
                    score = alignments[0].score
                else:
                    continue
                # below if statement determines what qualifies as a good alignment score
                if score >= key_n - key_n/10:
                    counter += 1
                    good_alignment_arr.append(alignments[0])

                    alignment_start = alignments[0].coordinates[0][0]
                    good_alignment_starting_pos.append(
                        str_dict[str][0]+alignment_start)

                    new_entry_start = str_dict[str][0]
                    new_entry_end = new_entry_start + alignment_start
                    new_str = self.input_s[new_entry_start:new_entry_end+1]
                    str_dict[new_str] = [new_entry_start, new_entry_end]
                    queue.append(new_str)

                    # new_entry_start = str_dict[str][0] + alignments[0].coordinates[0][0] + key_n + 1
                    new_entry_start = str_dict[str][0] + \
                        alignments[0].coordinates[0][-1] + 1
                    new_entry_end = str_dict[str][1]
                    new_str = self.input_s[new_entry_start:new_entry_end]
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

                    dict_insertions_and_deletions[str_dict[str][0]+alignment_start] = [
                        insertions_indexes, deletions_indexes]

                arr.append(score)
        return (good_alignment_starting_pos, arr, dict_insertions_and_deletions)

    def get_sorted_my_dict_keys(self):
        return sorted(self.my_dict, key=lambda k: self.my_dict[k].n_subStr)

    def filter_same_length(self):
        sorted_keys = self.get_sorted_my_dict_keys()
        prev_key = sorted_keys[0]
        for i in range(1, len(sorted_keys)):
            curr_key = sorted_keys[i]
            # if self.my_dict[prev_key].n_subStr == self.my_dict[curr_key].n_subStr and is_circular_rearrangement(prev_key, curr_key):
            if self.my_dict[prev_key].n_subStr == self.my_dict[curr_key].n_subStr:
                if len(self.my_dict[prev_key].indexes) + len(self.my_dict[prev_key].extra_alignment_indexes) >= len(self.my_dict[curr_key].indexes) + len(self.my_dict[prev_key].extra_alignment_indexes):
                    del self.my_dict[curr_key]
                else:
                    del self.my_dict[prev_key]
                    prev_key = curr_key
            else:
                prev_key = curr_key

    def get_if_n_subStr_equals_min_gap(self):
        output_arr = []
        for key in self.my_dict:
            if self.my_dict[key].n_subStr == self.my_dict[key].min_gap:
                output_arr.append(key)
        return output_arr

    def print_my_dict(self):
        # below is just for testing
        for key in self.my_dict:
            print(self.my_dict[key])
            print(key)
            print("\n")

    def print_my_dict_sorted(self):
        sorted_keys = self.get_sorted_my_dict_keys()

        for key in sorted_keys:
            print(self.my_dict[key])
            print(key)
            print("\n")

    def run(self):
        self.self_search()
        if len(self.my_dict.keys()) == 0:
            return 1

        self.expand_dict()
        if len(self.my_dict.keys()) == 0:
            return 1

        if not self.user_input_obj.ignore_alignment:
            for key in self.my_dict:
                self.my_dict[key].extra_alignment_indexes, self.my_dict[key].all_extra_alignment_scores, self.my_dict[
                    key].extra_alignment_insertions_and_deletions = self.alignment(key)
        self.filter_same_length()
        return 0


class full_analysis:

    def __init__(self, user_input_obj, all_chr_headers, same_prefix_length):
        self.user_input_obj = user_input_obj
        self.all_chr_headers = all_chr_headers
        self.same_prefix_length = same_prefix_length

    def graph_setup(self):
        self.fig, self.ax = plt.subplots(figsize=(16, 4))
        if self.user_input_obj.graph_title:
            self.ax.text(0, 9, self.user_input_obj.graph_title, ha='center',
                         va='center', fontsize=12, color='black')
        # ax.text(0, 8, "IT148", ha='center', va='center', fontsize=12, color='black')
        # ax.text(0, 7, "chr1", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 6.5, "chr2", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 6, "chr3", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 5.5, "chr4", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 5, "chr5", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 4.5, "chr6", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 4, "chr7", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 3.5, "chr8", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 3, "chr9", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 2.5, "chr10", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 2, "chr11", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 1.5, "chr12", ha='center', va='center', fontsize=10, color='black')
        # ax.text(0, 1, "chr13", ha='center', va='center', fontsize=10, color='black')

        # ax.text(0, 9, "IT148", ha='center', va='center', fontsize=12, color='black')
        self.ax.text(0, 8.5, "chr1", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 8, "chr2", ha='center', va='center',
                     fontsize=10, color='black')
        self.ax.text(0, 7.5, "chr3", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 7, "chr4", ha='center', va='center',
                     fontsize=10, color='black')
        self.ax.text(0, 6.5, "chr5", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 6, "chr6", ha='center', va='center',
                     fontsize=10, color='black')
        self.ax.text(0, 5.5, "chr7", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 5, "chr8", ha='center', va='center',
                     fontsize=10, color='black')
        self.ax.text(0, 4.5, "chr9", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 4, "chr10", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 3.5, "chr11", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 3, "chr12", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 2.5, "chr13", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 2, "chr14", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 1.5, "chr15", ha='center',
                     va='center', fontsize=10, color='black')
        self.ax.text(0, 1, "chr16", ha='center',
                     va='center', fontsize=10, color='black')

        # ax.set_xlabel("X-axis")
        # ax.set_ylim(.5, 8.5)
        self.ax.set_ylim(.5, 9.5)

        # Remove x/y-axis and add legend
        self.ax.get_yaxis().set_visible(False)
        self.ax.get_xaxis().set_visible(False)

    def graph_output(self, curr_sequence, chr_index):
        offset = 300
        counter = 0
        global_counter = chr_index
        max_x = 0
        min_x = 0
        # y_index = (14 - (global_counter//2)) / 2
        y_index = (17 - (global_counter//2)) / 2
        arrow_distance_orig = curr_sequence.n_subStr-1
        offset_offset = arrow_distance_orig//2
        set1 = copy.deepcopy(curr_sequence.indexes)
        set2 = copy.deepcopy(curr_sequence.extra_alignment_indexes)
        set2.sort()
        insertions_and_deletions = copy.deepcopy(
            curr_sequence.extra_alignment_insertions_and_deletions)
        all_points = []
        set1_ptr = 0
        set2_ptr = 0
        set1_n = len(set1)
        set2_n = len(set2)
        while set1_ptr < set1_n and set2_ptr < set2_n:
            set1_tmp = set1[set1_ptr]
            set2_tmp = set2[set2_ptr]
            if set1_tmp < set2_tmp:
                # all_points.append(set1_tmp + offset)
                all_points.append(set1_tmp + offset + (offset_offset*counter))
                set1_ptr += 1
            else:
                # all_points.append(set2_tmp + offset)
                all_points.append(set2_tmp + offset + (offset_offset*counter))
                set2_ptr += 1
            counter += 1
        while set2_ptr < set2_n:
            # all_points.append(set2[i] + offset)
            all_points.append(set2[set2_ptr] + offset +
                              (offset_offset * counter))
            set2_ptr += 1
            counter += 1
            # all_points.append(set1[i] + offset)
        while set1_ptr < set1_n:
            all_points.append(set1[set1_ptr] + offset +
                              (offset_offset*counter))
            set1_ptr += 1
            counter += 1

        # arrow_style = "->"
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
                insertions_and_deletions_key = abs(
                    all_points[i] - offset * sign)
                n_insertions = len(
                    insertions_and_deletions[insertions_and_deletions_key][0])
                n_deletions = len(
                    insertions_and_deletions[insertions_and_deletions_key][1])
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

            tail = Rectangle((bottom_left, y_index - tail_width / 2),
                             tail_length, tail_width, color=arrow_color)
            self.ax.add_patch(tail)

            if sign == -1 and point < min_x:
                min_x = point
            elif point > max_x:
                max_x = point

            triangle = Polygon([[point + (tail_length+head_length) * sign, y_index],
                                [point+(tail_length) * sign,
                                 y_index-head_width/2],
                                [point+(tail_length) * sign, y_index+head_width/2]],
                               closed=True, color=arrow_color)
            self.ax.add_patch(triangle)

            if not perfect_alignment:
                insertions_and_deletions_key = abs(
                    all_points[i] - offset * sign)
                insertions = insertions_and_deletions[insertions_and_deletions_key][0]
                deletions = insertions_and_deletions[insertions_and_deletions_key][1]
                for i in insertions:
                    # plt.plot([point+i*sign, point+i*sign], [y_index-.1, y_index+.1], color='blue', lw=1)  # Draw vertical line at midpoint
                    self.ax.plot([point+i*sign, point+i*sign], [y_index-.1,
                                 y_index+.1], color='blue', linestyle='-', lw=1)
                for i in deletions:
                    # plt.plot([point+i*sign, point+i*sign], [y_index-.1, y_index+.1], color='gold', lw=1)  # Draw vertical line at midpoint
                    self.ax.plot([point+i*sign, point+i*sign], [y_index-.1,
                                 y_index+.1], color='gold', linestyle='-', lw=1)

            # if previous_end_point is not None and abs(point - previous_end_point) > 1:
            if previous_end_point is not None and abs(point - previous_end_point) > offset_offset+1:
                # ax.plot([previous_end_point, point], [y_index, y_index], linewidth= 3, color='red', linestyle='-')
                self.ax.plot([previous_end_point+offset_offset/2*sign, point-offset_offset/2*sign], [
                             y_index, y_index], linewidth=2, color='red', linestyle='-')
            previous_end_point = end_point
            offset += offset_offset
        return max_x, min_x

    def save_graph(self, max_x, min_x):

        self.ax.set_xlim(min_x, max_x)
        self.fig.savefig(self.user_input_obj.graph_output,
                         dpi=self.user_input_obj.graph_dpi, bbox_inches='tight')

    def run_all_chr_ends(self, all_chr_ends):
        all_find_loops_objects = []
        no_loops_found_indexes = []
        counter = 0
        for i in range(len(all_chr_ends)):
            sequence = all_chr_ends[i]
            if not sequence:
                all_find_loops_objects.append(None)
                continue
            # print(f"same prefix length {i} = {self.same_prefix_length[i]}")
            loop_obj = find_loops(
                # Note: figure out why below line doesn't work correctly in all cases
                sequence[self.same_prefix_length[i]:], self.user_input_obj, self.same_prefix_length[i])
            # sequence, self.user_input_obj, self.same_prefix_length[i])
            err = loop_obj.run()
            # err code of 0 means that there were no issues
            if err == 0:
                all_find_loops_objects.append(loop_obj)
            else:
                no_loops_found_indexes.append(counter)
                all_find_loops_objects.append(None)
            counter += 1
        return all_find_loops_objects, no_loops_found_indexes

    def find_all_multi_chr_repeat_sequences(self, all_find_loops_objects):

        all_multi_chr_repeat_sequences_dict = defaultdict(list)
        for i in range(len(all_find_loops_objects)-1):
            # change curr name
            curr = all_find_loops_objects[i]
            if not curr:
                continue
            if self.user_input_obj.exaustive_analysis:
                possibilities = list(curr.my_dict.keys())
            else:
                possibilities = curr.get_if_n_subStr_equals_min_gap()
            possibilities.sort(key=len)
            for all_find_loops_objects_ptr in range(0, len(all_find_loops_objects)):
                if all_find_loops_objects_ptr == i or not all_find_loops_objects[all_find_loops_objects_ptr]:
                    continue
                # change curr_obj name
                curr_obj = all_find_loops_objects[all_find_loops_objects_ptr]
                sorted_my_dict_keys = curr_obj.get_sorted_my_dict_keys()
                possibilites_ptr = 0
                sorted_my_dict_keys_ptr = 0
                while possibilites_ptr < len(possibilities) and sorted_my_dict_keys_ptr < len(sorted_my_dict_keys):
                    if len(possibilities[possibilites_ptr]) > len(sorted_my_dict_keys[sorted_my_dict_keys_ptr]):
                        sorted_my_dict_keys_ptr += 1
                    elif len(possibilities[possibilites_ptr]) < len(sorted_my_dict_keys[sorted_my_dict_keys_ptr]):
                        possibilites_ptr += 1
                    elif is_circular_rearrangement(possibilities[possibilites_ptr], sorted_my_dict_keys[sorted_my_dict_keys_ptr]):
                        all_multi_chr_repeat_sequences_dict[(curr.my_dict[possibilities[possibilites_ptr]], i)].append(
                            (curr_obj.my_dict[sorted_my_dict_keys[sorted_my_dict_keys_ptr]], all_find_loops_objects_ptr))
                        del curr_obj.my_dict[sorted_my_dict_keys[sorted_my_dict_keys_ptr]]
                        # below += function is necessary
                        sorted_my_dict_keys_ptr += 1
                        possibilites_ptr += 1
                    else:
                        sorted_my_dict_keys_ptr += 1

        all_multi_chr_repeat_sequences_final_arr = []
        for key in all_multi_chr_repeat_sequences_dict:
            all_multi_chr_repeat_sequences_dict[key].append(key)
            all_multi_chr_repeat_sequences_dict[key].sort(key=lambda x: x[1])
            all_multi_chr_repeat_sequences_final_arr.append(
                all_multi_chr_repeat_sequences_dict[key])

        return all_multi_chr_repeat_sequences_final_arr

    def find_best_multi_chr_repeat_sequence(self, all_multi_chr_repeat_sequences_arr):
        most_chr_matches = 0
        best_chr_matches_index = 0
        for i in range(len(all_multi_chr_repeat_sequences_arr)):
            elem = all_multi_chr_repeat_sequences_arr[i]
            n = len(elem)
            if n > most_chr_matches:
                most_chr_matches = n
                best_chr_matches_index = i

        return all_multi_chr_repeat_sequences_arr[best_chr_matches_index]

    def run_full_analysis(self, all_chr_ends):
        all_find_loops_objects, no_loops_found_indexes = self.run_all_chr_ends(
            all_chr_ends)
        if len(all_find_loops_objects) == 0:
            print("no loops of the minimum size were found")
            return []
        all_multi_chr_repeat_sequences_final_arr = self.find_all_multi_chr_repeat_sequences(
            all_find_loops_objects)
        if len(all_multi_chr_repeat_sequences_final_arr) == 0:
            print("no loops were found in multiple chr's")
            return []
        best_repeat_sequence_arr = self.find_best_multi_chr_repeat_sequence(
            all_multi_chr_repeat_sequences_final_arr)
        if self.user_input_obj.graph_output:
            self.graph_setup()
            max_x = 0
            min_x = 0
            for elem in best_repeat_sequence_arr:
                max_x_output, min_x_output = self.graph_output(
                    elem[0], elem[1])
                max_x = max(max_x, max_x_output)
                min_x = min(min_x, min_x_output)
            self.save_graph(max_x, min_x)
        return best_repeat_sequence_arr, no_loops_found_indexes


class parse_fasta_file:

    def __init__(self, user_input_obj):
        self.user_input_obj = user_input_obj
        self.all_chr_headers = []
        self.all_chr_ends = []
        self.same_prefix_length = [
            user_input_obj.skip_prefix] * user_input_obj.maximum_ends

    # Function to extract the index from the header
    def get_index(self, header):
        # Use a regular expression to match digits at the end of the string
        match = re.search(r'(\d+)([LR])?$', header)
        if match:
            index = int(match.group(1))  # Extract the matched number
            # >>>>>>>>>>>>>>>>NOTe: see if this works when I don't have an
            # Extract 'L' or 'R', or default to empty
            modifier = match.group(2) if match.group(2) else ""
            modifier = match.group(2)
            if modifier:
                index *= 2
                if modifier == "L":
                    index -= 1

            if 1 <= index <= 32:
                return index - 1  # Convert to 0-based indexing
        raise ValueError(f"Invalid header format: {
                         header}, no sorting will be performed")

    def read_fasta_to_dict(self, file_path):
        fasta_dict = {record.id: str(record.seq)
                      for record in SeqIO.parse(file_path, "fasta")}
        return fasta_dict

    def read_original_reference_to_array(self):
        file_path = self.user_input_obj.original_reference
        sequences = []
        for sequence in SeqIO.parse(file_path, "fasta"):
            sequences.append(str(sequence.seq))
        return sequences

    def cmp_original_ref(self, original_ref_array):
        if len(original_ref_array) != self.user_input_obj.maximum_ends:
            print("number of entries in original reference file does not match the maximum ends value. Please change the maximum ends value with --maximum_ends flag or use a different reference file, origianl reference comparison will be ignored")

        for i in range(self.user_input_obj.maximum_ends):
            input_s = self.sequences_data[i]
            if not input_s:
                continue
            input_n = len(input_s)
            known_value = original_ref_array[i]
            known_value_n = len(known_value)

            output = 0

            for j in range(5):
                tmp_value = 0
                for k in range(min(known_value_n, input_n)):
                    if known_value[k] == input_s[k+j]:
                        tmp_value += 1
                    else:
                        break
                output = max(output, tmp_value)
                # if output > 5:
                # return(output+i)
            self.same_prefix_length[i] = max(
                self.same_prefix_length[i], output+j)

    def read_fasta_flexible(self):
        file_path = self.user_input_obj.file
        self.sequences_headers = [None] * self.user_input_obj.maximum_ends
        self.sequences_data = [None] * self.user_input_obj.maximum_ends
        total_fasta_entries_counter = 0
        with open(file_path, "r") as fasta_file:
            while True:
                counter = 0
                for Sequence in SeqIO.parse(fasta_file, "fasta"):
                    if counter >= self.user_input_obj.maximum_ends:
                        print("error, the number of entries in the fasta file exeeds the number of maximum ends. Please change the number of maximum ends withe the --maximum_ends flag")
                        return [], []
                    if self.user_input_obj.ignore_sorting:
                        self.sequences_headers[counter] = (str(Sequence.id))
                        self.sequences_data[counter] = (str(Sequence.seq))
                        counter += 1
                        total_fasta_entries_counter += 1
                    else:
                        try:
                            # Extract index from header
                            index = self.get_index(Sequence.id)
                            self.sequences_headers[index] = str(
                                Sequence.id)  # Assign sequence to the array
                            # Assign sequence to the array
                            self.sequences_data[index] = str(Sequence.seq)
                            total_fasta_entries_counter += 1
                        except ValueError as e:
                            print(e)
                            fasta_file.seek(0)
                            self.user_input_obj.ignore_sorting = True
                            total_fasta_entries_counter = 0
                            break
                else:
                    break

        is_whole_ref_sequence = False
        for elem in self.sequences_data:
            if elem and len(elem) > 100000:
                is_whole_ref_sequence = True
                break
        if is_whole_ref_sequence:
            print("whole reference detected\n")
        else:
            print("chr ends sequences detected\n")
        if is_whole_ref_sequence:
            if total_fasta_entries_counter > self.user_input_obj.maximum_ends//2:
                print("error, the number of entries in the fasta file exeeds the number of maximum ends for a whole reference sequence. Please change the number of maximum ends withe the --maximum_ends flag")
                return [], []
            # note: add an edge case for there being entry numbers that exceed the half of the maximum_ends variable
            for i in range(total_fasta_entries_counter-1, -1, -1):
                original_sequence = self.sequences_data[i]
                if original_sequence:
                    first_half_sequence = original_sequence[:5000]
                    second_half_sequence = original_sequence[-5000:]
                else:
                    first_half_sequence = None
                    second_half_sequence = None
                first_half_header = self.sequences_headers[i] + "_L"
                second_half_header = self.sequences_headers[i] + "_R"
                self.sequences_data[i*2+1] = second_half_sequence
                self.sequences_data[i*2] = first_half_sequence
                self.sequences_headers[i*2+1] = second_half_header
                self.sequences_headers[i*2] = first_half_header

    def find_flexible_telomeric_regions(self, window_size=50, threshold=0.9):
        # fix for full genome sequences
        all_chr_ends = self.sequences_data
        for i in range(len(all_chr_ends)):
            dna_sequence = all_chr_ends[i]
            if dna_sequence == None:
                continue
            regions = []
            end_regions = []
            n = len(dna_sequence)

            def is_ac_like(base1, base2):
                """Returns True if the pair is "AC-like" (e.g., AA, CC, CA, AC)."""
                return base1 in "AC" and base2 in "AC"

            def is_tg_like(base1, base2):
                """Returns True if the pair is "TG-like" (e.g., TT, GG, GT, TG)."""
                return base1 in "TG" and base2 in "TG"

            def is_dominated_by(window, motif_check):
                """Checks if a window is dominated by a given motif type."""
                valid_pairs = sum(1 for i in range(
                    len(window) - 1) if motif_check(window[i], window[i + 1]))
                return valid_pairs / (len(window) - 1) >= threshold

            def trim_to_valid_region(start, end, motif_type):
                # try to fix this if else statement to make it cleaner
                """Trims the region to ensure it starts and ends with valid characters."""
                # while start < end and dna_sequence[start] not in ("AC" if motif_type == "AC-like" else "TG"):
                # while start <= end-3 and dna_sequence[start] not in ("AC" if motif_type == "AC-like" else "TG"):
                # start += 1

                if motif_type == "AC":
                    while start <= len(dna_sequence) - 3:
                        if all(c in 'AC' for c in dna_sequence[start:start + 3]):
                            break
                        start += 1
                # while end > start and dna_sequence[end - 1] not in ("AC" if motif_type == "AC-like" else "TG"):
                    # end -= 1

                    while end >= start+2:
                        if all(c in 'AC' for c in dna_sequence[end - 3:end]):
                            break
                        end -= 1

                else:
                    while start <= len(dna_sequence) - 3:
                        if all(c in 'TG' for c in dna_sequence[start:start + 3]):
                            break
                        start += 1
                # while end > start and dna_sequence[end - 1] not in ("AC" if motif_type == "AC-like" else "TG"):
                    # end -= 1

                    while end >= start+2:
                        if all(c in 'TG' for c in dna_sequence[end - 3:end]):
                            break
                        end -= 1

                return start, end

            for j in range(n - window_size + 1):
                window = dna_sequence[j:j + window_size]
                if is_dominated_by(window, is_ac_like):
                    regions.append((j, j + window_size, "AC"))
                elif is_dominated_by(window, is_tg_like):
                    regions.append((j, j + window_size, "TG"))

            # Merge and trim overlapping regions of the same motif type
            merged_regions = []
            for start, end, motif_type in regions:
                if not merged_regions or start > merged_regions[-1][1] or motif_type != merged_regions[-1][2]:
                    merged_regions.append((start, end, motif_type))
                else:
                    merged_regions[-1] = (merged_regions[-1]
                                          [0], end, motif_type)

            # Trim regions to valid start and end
            trimmed_regions = [
                (start, end, motif_type)
                for start, end, motif_type in (trim_to_valid_region(s, e, t) + (t,) for s, e, t in merged_regions)
                if start < end  # Ensure trimmed region is valid
            ]

            valid_telomer = None
            for elem in trimmed_regions:
                if elem[1] - elem[0] >= 100 and (elem[0] <= 200 or elem[1] >= n-200):
                    if valid_telomer:
                        print(
                            f"multiple valid sequences were found in index {i}")
                        if elem[1] - elem[0] > len(valid_telomer):
                            valid_telomer = dna_sequence[elem[0]:elem[1] + 1]
                            if elem[2] == "AC":
                                valid_telomer = str(
                                    Seq(valid_telomer).reverse_complement())
                    else:
                        valid_telomer = dna_sequence[elem[0]:elem[1] + 1]
                        if elem[2] == "AC":
                            valid_telomer = str(
                                Seq(valid_telomer).reverse_complement())

            all_chr_ends[i] = valid_telomer
        return all_chr_ends

    def run(self):
        # all_chr_ends = self.read_fasta_to_array(file_path)
        self.read_fasta_flexible()
        self.find_flexible_telomeric_regions()
        if self.user_input_obj.original_reference:
            original_ref_arr = self.read_original_reference_to_array()
            self.cmp_original_ref(original_ref_arr)


def mod_str(my_str, self_search_type_object):
    # make this function less shit
    offset = self_search_type_object.n_subStr
    perfect_indexes_arr = self_search_type_object.indexes
    extra_alignment_indexes_arr = self_search_type_object.extra_alignment_indexes
    extra_alignment_indexes_arr.sort()
    # extra_alignment_insertions_and_deletions_dict = self_search_type_object.extra_alignment_insertions_and_deletions
    perfect_indexes_arr_ptr = len(perfect_indexes_arr)-1
    extra_alignment_indexes_arr_ptr = len(extra_alignment_indexes_arr)-1
    # for i in range(len(my_dict)-1, -1, -1):
    # my_str = my_str[:my_dict[i]] + "[" + my_str[my_dict[i]:my_dict[i] + offset] + "]" + my_str[my_dict[i]+offset:]
    while perfect_indexes_arr_ptr >= 0 and extra_alignment_indexes_arr_ptr >= 0:
        if perfect_indexes_arr[perfect_indexes_arr_ptr] >= extra_alignment_indexes_arr[extra_alignment_indexes_arr_ptr]:
            curr_index = perfect_indexes_arr[perfect_indexes_arr_ptr]
            perfect_indexes_arr_ptr -= 1
            perfect_alignment = True
            right_delimeter = "]"
            left_delimeter = "["
        else:
            curr_index = extra_alignment_indexes_arr[extra_alignment_indexes_arr_ptr]
            extra_alignment_indexes_arr_ptr -= 1
            perfect_alignment = False
            right_delimeter = "}"
            left_delimeter = "{"
        my_str = my_str[:curr_index+offset] + \
            right_delimeter + my_str[curr_index+offset:]
        if not perfect_alignment:
            insertions = self_search_type_object.extra_alignment_insertions_and_deletions[
                curr_index - self_search_type_object.offset][0]
            deletions = self_search_type_object.extra_alignment_insertions_and_deletions[
                curr_index - self_search_type_object.offset][1]
            insertions_ptr = 0
            deletions_ptr = 0
            n_insertions = len(insertions)
            n_deletions = len(deletions)
            while insertions_ptr < n_insertions and deletions_ptr < n_deletions:
                if insertions[insertions_ptr] < deletions[deletions_ptr]:
                    my_str = my_str[:curr_index + insertions[insertions_ptr]] + \
                        "^" + my_str[curr_index + insertions[insertions_ptr]:]
                    insertions_ptr += 1
                else:
                    my_str = my_str[:curr_index + deletions[deletions_ptr]] + \
                        "!" + my_str[curr_index + deletions[deletions_ptr]:]
                    deletions_ptr += 1
            while insertions_ptr < n_insertions:
                my_str = my_str[:curr_index + insertions[insertions_ptr]] + \
                    "^" + my_str[curr_index + insertions[insertions_ptr]:]
                insertions_ptr += 1
            while deletions_ptr < n_deletions:
                my_str = my_str[:curr_index + deletions[deletions_ptr]] + \
                    "!" + my_str[curr_index + deletions[deletions_ptr]:]
                deletions_ptr += 1
        my_str = my_str[:curr_index] + left_delimeter + my_str[curr_index:]
    while perfect_indexes_arr_ptr >= 0:
        curr_index = perfect_indexes_arr[perfect_indexes_arr_ptr]
        my_str = my_str[:curr_index+offset] + "]" + my_str[curr_index+offset:]
        my_str = my_str[:curr_index] + "[" + my_str[curr_index:]
        perfect_indexes_arr_ptr -= 1
    while extra_alignment_indexes_arr_ptr >= 0:
        curr_index = extra_alignment_indexes_arr[extra_alignment_indexes_arr_ptr]
        my_str = my_str[:curr_index+offset] + "}" + my_str[curr_index+offset:]
        insertions = self_search_type_object.extra_alignment_insertions_and_deletions[
            curr_index - self_search_type_object.offset][0]
        deletions = self_search_type_object.extra_alignment_insertions_and_deletions[
            curr_index - self_search_type_object.offset][1]
        insertions_ptr = 0
        deletions_ptr = 0
        n_insertions = len(insertions)
        n_deletions = len(deletions)
        while insertions_ptr < n_insertions and deletions_ptr < n_deletions:
            if insertions[insertions_ptr] < deletions[deletions_ptr]:
                my_str = my_str[:curr_index + insertions[insertions_ptr]] + \
                    "^" + my_str[curr_index + insertions[insertions_ptr]:]
                insertions_ptr += 1
            else:
                my_str = my_str[:curr_index + deletions[deletions_ptr]] + \
                    "!" + my_str[curr_index + deletions[deletions_ptr]:]
                deletions_ptr += 1
        while insertions_ptr < n_insertions:
            my_str = my_str[:curr_index + insertions[insertions_ptr]] + \
                "^" + my_str[curr_index + insertions[insertions_ptr]:]
            insertions_ptr += 1
        while deletions_ptr < n_deletions:
            my_str = my_str[:curr_index + deletions[deletions_ptr]] + \
                "!" + my_str[curr_index + deletions[deletions_ptr]:]
            deletions_ptr += 1
        my_str = my_str[:curr_index] + "{" + my_str[curr_index:]
        extra_alignment_indexes_arr_ptr -= 1

    my_str = my_str[:self_search_type_object.offset] + \
        "|" + my_str[self_search_type_object.offset:]

    return my_str


def main(args):
    # converts args to a dictionary and pass it to the user_input class
    user_input_obj = user_input(**vars(args))
    # Check if the file exists
    if not os.path.isfile(user_input_obj.file):
        print(f"Error: The file '{user_input_obj.file}' does not exist.")
        return
    parse_fasta_file_obj = parse_fasta_file(user_input_obj)
    # all_chr_headers, all_chr_ends = parse_fasta_file_obj.run("AAAAAallFiles_145_modified.txt")
    parse_fasta_file_obj.run()
    all_chr_headers = parse_fasta_file_obj.sequences_headers
    all_chr_ends = parse_fasta_file_obj.sequences_data
    same_prefix_length = parse_fasta_file_obj.same_prefix_length
    if len(all_chr_ends) == 0:
        return
    full_analysis_obj = full_analysis(
        user_input_obj, all_chr_headers, same_prefix_length)
    best_repeat_sequence_arr, no_loops_found_indexes = full_analysis_obj.run_full_analysis(
        all_chr_ends)
    # for elem in best_repeat_sequence_arr:
    best_repeat_sequence_arr_ptr = 0
    no_loops_found_indexes_ptr = 0
    n_no_loops_found_indexes = len(no_loops_found_indexes)
    n_best_repeat_sequences_arr = len(best_repeat_sequence_arr)
    elem = best_repeat_sequence_arr[0]
    sequence = all_chr_ends[elem[1]][elem[0].indexes[0] +
                                     elem[0].offset:elem[0].indexes[0] + elem[0].n_subStr + elem[0].offset]
    print(f"best repeating sequence: {sequence}\n")
    printed_sequnce = False
    for i in range(len(all_chr_headers)):
        # if elem[1] == i:
        if best_repeat_sequence_arr_ptr < n_best_repeat_sequences_arr and best_repeat_sequence_arr[best_repeat_sequence_arr_ptr][1] == i:
            elem = best_repeat_sequence_arr[best_repeat_sequence_arr_ptr]
            print(all_chr_headers[elem[1]])
            # print(f"chr{convert_num_to_chr_end(elem[1])}")
            elem[0].add_offset()
            print(elem[0])
            moded_str = mod_str(all_chr_ends[elem[1]], elem[0])
            print(f"{moded_str}\n")
            best_repeat_sequence_arr_ptr += 1
        elif no_loops_found_indexes_ptr < n_no_loops_found_indexes and i == no_loops_found_indexes[no_loops_found_indexes_ptr]:
            print(f"{all_chr_headers[i]}\nno loops of the minimum length were found in {
                  all_chr_headers[i]}\n")
            no_loops_found_indexes_ptr += 1
        else:
            print(f"{all_chr_headers[i]}\nthe repeated sequence was not found in {
                  all_chr_headers[i]}\n")
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="program used for finding repeating sequences inside of data, provided in the form of a .fasta file")
    parser.add_argument("file", help="The path to the file to be processed.")
    parser.add_argument(
        "-o", "--output", help="Optional output file to save the content.")
    parser.add_argument("--min_length", type=int, default=50,
                        help="The minimum length for a valid repeat sequence (default: 50)")
    parser.add_argument("--graph_dpi", type=int, default=300,
                        help="the dpi of the saved graph (default: 300)")
    parser.add_argument("-go", "--graph_output",
                        help="Optional output file that a graph of the output will be saved too. (if No output is given, no graph will be created)")
    parser.add_argument("-ea", "--exaustive_analysis", action="store_true",
                        help="performs a more exaustive analysis which sacrifices runtime for the chance to catch some edgecase sequences")
    parser.add_argument("-me", "--maximum_ends", type=int, default=32,
                        help="changes the maximum number of chr ends that are expected (default is 32)")
    parser.add_argument("-is", "--ignore_sorting", action="store_true",
                        help="stops the program not attempt to automatically sort the sequences from fasta file based on the sequence headers")
    parser.add_argument("-ia", "--ignore_alignment", action="store_true",
                        help="stops the program from attempting to find imperfecting alighments. This drastically decreases the runtime of the program")
    parser.add_argument("-or", "--original_reference",
                        help="Optional input file with the original reference. New chr ends will be compared against the original reference, and any identical prefixes will not be examined for repeat sequences")
    parser.add_argument("-sp", "--skip_prefix", type=int, default=0,
                        help="skips looking for repeating sequences for the specified number of base pairs in all chr's")
    parser.add_argument(
        "-gt", "--graph_title", help="optional input for graph title (must be used with --graph_output flag)")
    args = parser.parse_args()

    # Open the output file if provided
    with open(args.output, 'w') if args.output else sys.stdout as output:
        # Redirect all prints within main
        with redirect_stdout(output):
            main(args)
