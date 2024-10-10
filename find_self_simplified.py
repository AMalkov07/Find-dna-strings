from collections import deque

class self_search_type:
    def __init__(self, n_subStr, indexes):
        self.n_subStr = n_subStr
        #self.subStr_start = subStr_start
        self.indexes = indexes
        #self.coverage = [indexes[0], indexes[0]+n_subStr]

    def __repr__(self):
        n = self.n_subStr
        return f"substring length={self.n_subStr}, number of matches={len(self.indexes)}\n{self.indexes}\n"

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
                #print(f"start: {start}")
                #print(f"end: {end}")
                #print(f"cmp_start: {cmp_start}")
                #print(f"cmp_end: {cmp_end}")
                del my_dict[input_s[cmp_start:cmp_end+1]]
                #print(f"{input_s[cmp_start:cmp_end+1]} has been deleted")
                #if end > cmp_end:
                    #del my_dict[input_s[cmp_start:cmp_end+1]]
                #else:
                    #del my_dict[input_s[start:end+1]]
    print(my_dict)

                


def expand_dict(input_s, my_dict, coverage_dict):
    input_s_n = len(input_s)
    queue = deque()
    for key in my_dict.keys():
        queue.append(key)
    while queue:
        #print(f"my_dict: {my_dict}")
        #print(f"queue: {queue}")
        #print(f"coverage_dict: {coverage_dict}")
        curr = queue.popleft()
        n_string = my_dict[curr].n_subStr
        indexes = my_dict[curr].indexes
        n_indexes = len(indexes)
        #counter = -1
        s = ""
        tmp_dict = {}
        for index in indexes:
            r_index = index + n_string
            if r_index >= input_s_n:
                continue
            s = input_s[index:r_index+1]
            if s in my_dict:
                my_dict[s].indexes.append(index)
                tmp_dict[s] = index
                #counter += 1
            elif s in tmp_dict:
                my_dict[s] = self_search_type(len(s), [tmp_dict[s], index])
                queue.append(s)
                #counter = 2
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
                        del my_dict[input_s[start_index:end_index+1]]
                    coverage_dict[n_tmp_indexes][start_index] = start_index + tmp.n_subStr - 1
    consolidate(input_s, my_dict, coverage_dict)

            


def self_search(input_s):
    my_dict = {}
    coverage_dict = {}
    n = len(input_s)
    min_length = 50
    #min_length = 10
    # i controlls the starting location
    for i in range(0, len(input_s)-min_length+1):
        #j controlls the size of the search window
        #for j in range(min_length, n-i):
        j = min_length
        subStr = input_s[i:i+j]
        if subStr in my_dict:
            continue
        # k controlls the sliding
        # using i+1 means that we will find strings that are offsets of the same string
        counter = 1
        for k in range(i+1, n-j+1):
            if input_s[k:k+j] == subStr:
                counter += 1
                if subStr in my_dict:
                    my_dict[subStr].indexes.append(k)
                else:
                    my_dict[subStr] = self_search_type(j, [i,k])
        if counter > 1:
            # coverage_dict is a dictionary of dictionaries. First value refers to how many matches for a particular value, second value is the start of the coverage, and the dictionary value is the end of the coverage
            if counter not in coverage_dict:
                coverage_dict[counter] = {}
            coverage_dict[counter][i] = i+j - 1

    expand_dict(input_s, my_dict, coverage_dict)

#self_search("abcdwabcdwbcdwabcd")


#48 - chr1l
#self_search("GTGTGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGGGTGGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGGTGTGGTGTGTGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGTGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGGGTGTGGGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGGGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGTGTGGGGTGTGGTGTGTGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGGGGTGTGGTGTGTGTGTGTGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTG")

#48 - chr1r
#self_search("GGTGTGTGTGGGTGTGGTGTGGGTGTGGTTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGTGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGTGGGTGTGGTGTGTGTGTGGTGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGGGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGTGGGTGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGTGTGTGTGTGTGTGTGTGGGTGTGTGTGTGTGTGGTGTGTGTGTGTGGGTGTGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGTGTGGGTGTGTGGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGGTGTGGTGTGTGTGTGGGTGTGGGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTG")

#48 - chr2l
#self_search("ACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACCCACACACCACACCCACACACCCACACACCACACCCACACACCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACCCCACACCCACACACCACACCCACACACCCACACACCCACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCACACCCACACACACACCACACCCCACACACACACCACACCCACACACCCACACCCACACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACACCACACACACCACACCCACACACACACACACCCACACACCCACACACCACACCCACACACACACACACACCACACACACACACACACACACACCACACACCCACACACCACACCCACACACACCCACACCACACACACACACACCACACACACACACACCCACACACCACACCCACACACCACACCCACACACCCACACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACACCACACACACCACACCCACACACCCACACCCACACCCACACACACACCCCACACACCACACCCACACACCCACACCCACACCACACCCACACACACCACACACCACACCCACACACCACACCCACACACCACCCACACCCACACAC")

#48 - chr2r
self_search("ACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACCCACACACCACACCCACACACCCACACACCACACCCACACACCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACCCCACACCCACACACCACACCCACACACCCACACACCCACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCACACCCACACACACACCACACCCCACACACACACCACACCCACACACCCACACCCACACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCACACCCACACACCACACACACCACACCCACACACACACACACCCACACACCCACACACCACACCCACACACACACACACACCACACACACACACACACACACACCACACACCCACACACCACACCCACACACACCCACACCACACACACACACACCACACACACACACACCCACACACCACACCCACACACCACACCCACACACCCACACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACACCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACCCACACACACCACACCCACACACCCACACCCACACCCACACACACACCACACCCCACACACACACCACACCCACACCACACCCACACACCCACACACCACACCCACACACCACACCCACACACCCACACACCACACACCACACACACCACACCCACACACCCACACCCACACCCACACACACACCCCACACACCACACCCACACACCCACACCCACACCACACCCACACACACCACACACCACACCCACACACCACACCCACACACCACCCACACCCACACAC")