class self_search_type:
    def __init__(self, n_subStr, matches_start_indexes):
        self.n_subStr = n_subStr
        #self.subStr_start = subStr_start
        self.matches_start_indexes = matches_start_indexes

    def __repr__(self):
        return f"substring length={self.n_subStr}, number of matches={len(self.matches_start_indexes)}\n"

def expand_dict(my_dict):
    for key in my_dict.keys():
        print("hi")

def self_search(input_s):
    my_dict = {}
    n = len(input_s)
    min_length = 2 
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
        for k in range(i+1, n-j+1):
            if input_s[k:k+j] == subStr:
                if subStr in my_dict:
                    my_dict[subStr].matches_start_indexes.append(k)
                else:
                    my_dict[subStr] = self_search_type(j, [i,k])

    print(my_dict)
    #sorted_dict = sorted(my_dict.items(), key=lambda item: item[1].n_subStr, reverse=True)
    #print(sorted_dict)

self_search("abcdabc")