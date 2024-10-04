####
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
        # Add last entry
        if identifier:
            fasta_dict[identifier] = ''.join(sequence_lines)
    return fasta_dict

#file_path = r"C:\Users\Andrey\Desktop\momWork\2024\Find-dna-strings\6991_only_telomeres_new.fasta"
file_path = r"6991_only_telomeres_new.fasta"

dict = read_fasta(file_path)

#class used for the creation of a custom variable type that will store all the info that we are interested in
class customType:
    def __init__(self, input_s_li, input_s_ri, dict_str, key_index, is_reverse):
        self.input_s_li = input_s_li
        self.input_s_ri = input_s_ri
        self.dict_str = dict_str
        self.key_index = key_index
        self.is_reverse = is_reverse
    
    def __repr__(self):
        return f"CustomType(input_s_li={self.input_s_li}, input_s_ri={self.input_s_ri}, dict_str={self.dict_str}, key_index={self.key_index}, is_reverse={self.is_reverse})"
    
# pairwise comparison

from Bio.Align import PairwiseAligner

# create the aligner
aligner = PairwiseAligner()

# Customize the alignment parameters
aligner.mode = 'local'  # Local (Smith-Waterman) alignment; use 'global' for global (Needleman-Wunsch)
aligner.match_score = 2  # Score for a match
aligner.mismatch_score = -1  # Penalty for a mismatch
#aligner.open_gap_score = -1  # Penalty for opening a gap
aligner.open_gap_score = -1  # Penalty for opening a gap
aligner.extend_gap_score = -0.1  # Penalty for extending a gap

input_s = "GTGTGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGGTGTGTGTGTGGGTGTGTGGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGTGGGTGTGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGGGTGTGGTGTGTGGGTGTGG"

max_length = 250
output = {}

# main function that finds whether or not parts of out input_s string are contained inside of our dictionary strings as substrings
def rec(input_s, offset):
    #print(f"running rec on {input_s}")
    if len(input_s) < 20:
        return
    for i in range(max_length, -1, -1):
        sub_str = input_s[0:i]
        for key in dict.keys():
            alignments = aligner.align(sub_str, dict[key])
            #for alignment in alignments:
            alignment = alignments[0]
            if alignment.score >= i*2 - 2:
                print(key)
                print(i)
                print(f"Alignment score: {alignment.score}")
                print(alignment)
                print(f"start: {offset}")
                print(f"end: {offset+i}")
                print("\n")
                rec(input_s[i:], offset+i)
                return
    return
            #if dict[key]== sub_str:
            #if sub_str in dict[key]: 
                #index = dict[key].find(sub_str)
                #output[sub_str] = customType(offset, offset+i, key, index, False)
                #rec(input_s[i:], offset+i)
                #return
            


class self_search_type:
    def __init__(self, n_subStr, subStr_start, matches_start_indexes):
        self.n_subStr = n_subStr
        self.subStr_start = subStr_start
        self.matches_start_indexes = matches_start_indexes

    def __repr__(self):
        return f"substring length={self.n_subStr}, number of matches={len(self.matches_start_indexes)}\n"

def self_search(input_s):
    my_dict = {}
    n = len(input_s)
    min_length = 10 
    for i in range(0, len(input_s)//2):
        for j in range(min_length, n-i):
            subStr = input_s[i:i+j]
            if subStr in my_dict:
                continue
            for k in range(i+1, n-j+1):
                if input_s[k:k+j] == subStr:
                    if subStr in my_dict:
                        my_dict[subStr].matches_start_indexes.append(k)
                    else:
                        my_dict[subStr] = self_search_type(j, i, [k])
    sorted_dict = sorted(my_dict.items(), key=lambda item: item[1].n_subStr, reverse=True)
    print(sorted_dict)

#rec(input_s, 0)

def compare_start(chr_compare, input_s):
    compare_str = dict[chr_compare]
    n = len(compare_str)
    counter = 10
    alignments = aligner.align(input_s[0:counter], compare_str[0:counter])
    alignment = alignments[0]
    # max mistakes currently coded as 2
    while alignment.score >= counter * 2 -(2):
    #while alignment.score >= counter * 2 -(1+counter//10):
        counter += 1
        alignments = aligner.align(input_s[0:counter], compare_str[0:counter])
        alignment = alignments[0]

    print(counter)
    print(alignment)
    print(alignment.score)

    print(">>>>>>>")
    print(compare_str)
    print(input_s)


#compare_start("chr1L_Telomere_Repeat", input_s)
#self_search(input_s)


compare_str = dict["chr1L_Telomere_Repeat"]
n = len(compare_str)
alignments = aligner.align(input_s[0:n+5], compare_str)
#for alignment in alignments:
#for i in range(len(alignments)):
print(len(alignments))
alignment = alignments[1000]
print(f"score: {alignment.score}")
print(alignment)
print(">>>>>")