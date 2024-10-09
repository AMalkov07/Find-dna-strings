from collections import deque, defaultdict

class self_search_type:
    def __init__(self, n_subStr, indexes):
        self.n_subStr = n_subStr
        #self.subStr_start = subStr_start
        self.indexes = indexes
        self.coverage = [indexes[0], indexes[0]+n_subStr]

    def __repr__(self):
        return f"substring length={self.n_subStr}, number of matches={len(self.indexes)}\n"

def expand_dict(input_s, my_dict, coverage_dict):
    input_s_n = len(input_s)
    queue = deque()
    for key in my_dict.keys():
        queue.append(key)
    while queue:
        print(queue)
        print(coverage_dict)
        tmp = queue.popleft()
        n_string = my_dict[tmp].n_subStr
        indexes = my_dict[tmp].indexes
        n_indexes = len(indexes)
        if indexes[0] not in coverage_dict[len(indexes)]:
            del my_dict[tmp]
            continue
        counter = -1
        s = ""
        tmp_dict = {}
        for index in indexes:
            r_index = index + n_string
            if r_index >= input_s_n:
                continue
            s = input_s[index:r_index+1]
            if s in my_dict:
                my_dict[s].indexes.append(index)
                counter += 1
            if s in tmp_dict:
                my_dict[s] = self_search_type(len(s), [tmp_dict[s], index])
                queue.append(s)
                counter = 2
            else:
                tmp_dict[s] = index
        if counter > 0:
            r_index = indexes[0] + n_string
            if counter not in coverage_dict:
                coverage_dict[counter] = {indexes[0]: r_index}
            elif r_index in coverage_dict[counter]:
                coverage_dict[counter][indexes[0]] = coverage_dict[counter][r_index]
                del coverage_dict[counter][r_index]
            else:
                coverage_dict[counter][indexes[0]] = r_index

        if counter == n_indexes:
            #print(my_dict)
            #print(s)
            #print(len(s))
            #print(counter)
            del my_dict[s[:-1]]
    print(my_dict)





def self_search(input_s):
    #my_dict = defaultdict(lambda: 0)
    my_dict = {}
    coverage_dict = {}
    n = len(input_s)
    min_length = 2 
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
    #sorted_dict = sorted(my_dict.items(), key=lambda item: item[1].n_subStr, reverse=True)
    #print(sorted_dict)

#self_search("GTGTGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGG")
self_search("abcdwbcwabcd")