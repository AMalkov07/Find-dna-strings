from collections import deque, defaultdict
import sys
from Bio import Align
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon
import copy

def is_circular_rearrangement(s1, s2):
    # Check if the strings are of the same length
    if len(s1) != len(s2):
        return False
    
    # Concatenate s2 with itself
    s2_double = s2 + s2
    
    # Check if s1 is a substring of s2_double
    return s1 in s2_double

class self_search_type:
    def __init__(self, n_subStr, indexes, min_gap, has_overlap):
        self.n_subStr = n_subStr
        self.indexes = indexes
        self.min_gap = min_gap
        self.has_overlap = has_overlap
        self.n_extra_alignment = 0
        self.extra_alignment_indexes = []
        self.extra_alignment_insertions_and_deletions = {}
    
    #def average_distance(self):
        #arr = self.indexes
        #tmp = 0
        #n = len(arr)
        #for i in range(1, n):
            #tmp += (arr[i] - arr[i-1])
        #return tmp//n

    
    def __repr__(self):
        n = self.n_subStr
        #avg_dist = self.average_distance()
        return f"substring length={self.n_subStr}, number of matches={len(self.indexes)}, min_gap={self.min_gap}, overlap={self.has_overlap}, number of extra alignments={self.n_extra_alignment}\nstarting indexes: {self.indexes}\nextra alignemtn indexes: {self.extra_alignment_indexes}\n"

class find_loops:

    def __init__(self, input_s):
        self.input_s = input_s
        self.my_dict = {}
        self.coverage_dict = defaultdict(dict)

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
            #n_indexes = len(indexes)
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
                    my_dict[s] = self_search_type(len(s), [first_val, index], gap, has_overlap)
                    queue.append(s)
                else:
                    tmp_dict[s] = index
            for key in tmp_dict.keys():
                #if key in my_dict:
                    #tmp = my_dict[key]
                    #n_tmp_indexes = len(tmp.indexes)
                    #start_index = tmp.indexes[0]
                    #if n_tmp_indexes not in coverage_dict:
                        #coverage_dict[n_tmp_indexes] = {start_index: start_index + tmp.n_subStr - 1}
                    #else:
                        #if  start_index in coverage_dict[n_tmp_indexes]:
                            #end_index = coverage_dict[n_tmp_indexes][start_index]
                            #sub_str = input_s[start_index:end_index+1]
                            #if my_dict[sub_str].is_overlap == my_dict[key].is_overlap:
                                #del my_dict[input_s[start_index:end_index+1]]
                        #coverage_dict[n_tmp_indexes][start_index] = start_index + tmp.n_subStr - 1
                if key in my_dict:
                    # if statement used for imediately removing sequences with overlap (added by previous loop)
                    if not my_dict[key].has_overlap:
                        curr_dict_value = my_dict[key]
                        n_curr_dict_value_indexes = len(curr_dict_value.indexes) # remove this variable
                        if n_curr_dict_value_indexes == len(my_dict[curr].indexes):
                            del my_dict[curr]

    def self_search(self, input_s = "", min_length = 50):
        #self.input_s = input_s
        input_s = self.input_s
        n = len(input_s)
        tmp_dict = defaultdict(list)

        # this loops creates tmp_dict which contains all possible substrings of size min_length as the dictionary key, and the indexes of that substring as the value
        for i in range(0, n-min_length+1):
            sub_str = input_s[i:i+min_length]
            tmp_dict[sub_str].append(i)

        # coverage_dict is a dictionary of dictionaries. First value refers to how many matches for a particular value, second value is the start of the coverage, and the dictionary value is the end of the coverage

        # following code is for converting 
        last_index = None
        last_key = None
        should_skip = False
        sequential_englargments_counter = 1

        for i in range(0, n-min_length+1):
            sub_str = input_s[i:i+min_length]
            tmp_arr = tmp_dict[sub_str]
            len_arr = len(tmp_arr)

            # if the array from our tmp_dict value only contains 1 value then that means it didn't repeat at all and we are not interested in it.
            # following if conditional also determines the minimum number of matches that must be found in order to continue
            if len_arr > 2:

                # following if statement is for expading the substrings that we originally found from the min_length to a larger length, if we can easily tell that the found repeating string should be longer the min_length
                if last_index and last_index == i-1:
                    prev_arr = self.my_dict[last_key].indexes
                    len_prev = len(prev_arr)
                    has_all_same_locations = True

                    if len_arr == len_prev:
                        for j in range(0,len_prev):
                            if tmp_arr[j] != prev_arr[j] + sequential_englargments_counter:
                                has_all_same_locations = False
                                break

                        if has_all_same_locations:
                            should_skip = True
                            #tmp_object = self.my_dict[last_key]
                            tmp_object = copy.deepcopy(self.my_dict[last_key])
                            had_overlapping = tmp_object.min_gap < tmp_object.n_subStr
                            tmp_object.n_subStr += 1
                            tmp_object.has_overlap = tmp_object.min_gap < tmp_object.n_subStr
                            new_key = str(last_key)+sub_str[-1]
                            self.my_dict[new_key] = tmp_object
                            if had_overlapping == tmp_object.has_overlap:
                                del self.my_dict[last_key]
                            last_index = i
                            last_key = new_key
                            sequential_englargments_counter += 1

                # following if statement is for converting the basic array that serves as the value to our tmp_dict dictionary, to a more robust object type that can tonain more information, and we store this object type in the my_dict dictionary
                if not should_skip:
                    min_gap = sys.maxsize
                    for j in range(1, len_arr):
                        min_gap = min(min_gap, tmp_arr[j] - tmp_arr[j-1])
                    first_index = tmp_arr[0]
                    self.coverage_dict[len_arr][first_index] = first_index+min_length - 1
                    overlap = min_gap < min_length
                    self.my_dict[sub_str] = self_search_type(min_length, tmp_arr, min_gap, overlap)
                    last_index = i
                    last_key = sub_str
                    sequential_englargments_counter = 1
                should_skip = False

        if len(list(self.my_dict.keys())) == 0:
            print("no loops of the minimum length were found\n")
            return


    def print_my_dict(self):
        # below is just for testing
        for key in self.my_dict:
            print(self.my_dict[key])
            print(key)
            print("\n")

    def run(self):
        self.self_search()
        self.expand_dict()
        self.print_my_dict()

def main():

    my_chr_end = find_loops("ACCCACACCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACCACACCCACACACACACACCCACACACACACACCCACACCACACCCACACACACACCCACACCACACCCACACCACACCACACCCACACACACACCCACACACACACCCACACCACACCCACACCCACACACCACACCCACACCCACACACACCACACCCACACCCACACCCACACACCACACCCACACCCACACACCACACCCACACCCACACCCACACACCCACACCACACCCACACACCACACCCACACCACACCCACACACCCACACACACACACCCACACACCACACCCACACCACACCCACACCACACCCACACCCACACCCACACCACACCCACACACACCACACCCACACCCACACACCACACACAC")
    my_chr_end.run()


if __name__ == "__main__":
    main()
