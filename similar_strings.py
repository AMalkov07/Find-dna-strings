#dict = {"s1": "cat",
     #"s2": "dog",
     #"s3": "cow"}


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

file_path = r"C:\Users\Andrey\Desktop\momWork\2024\Find-dna-strings\6991_only_telomeres_new.fasta"

dict = read_fasta(file_path)

class customType:
    def __init__(self, input_s_li, input_s_ri, dict_str, key_index, is_reverse):
        self.input_s_li = input_s_li
        self.input_s_ri = input_s_ri
        self.dict_str = dict_str
        self.key_index = key_index
        self.is_reverse = is_reverse
    
    def __repr__(self):
        return f"CustomType(input_s_li={self.input_s_li}, input_s_ri={self.input_s_ri}, dict_str={self.dict_str}, key_index={self.key_index}, is_reverse={self.is_reverse})"
    
input_s = "TGTGTGTGGTGTGTGGGTGTGGGTGTGGTGTGTGTGGGTGGGTGTGGTGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGTGTGTGTGGGTGTGTGGGGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGGGTGTGGTGTGTGTGGGTGTGGTGTGTGGGTGTGGTGTGTGGGTGTGGGTGTGTGGGTGTGGGTG"

max_length = 60
output = {}

def rec(input_s, offset):
    #print(f"running rec on {input_s}")
    if len(input_s) < 20:
        return
    for i in range(max_length, -1, -1):
        sub_str =input_s[0:i]
        for key in dict.keys():
            #if dict[key]== sub_str:
            if sub_str in dict[key]: 
                index = dict[key].find(sub_str)
                output[sub_str] = customType(offset, offset+i, key, index, False)
                rec(input_s[i:], offset+i)
                return


rec(input_s, 0)
for key in output.keys():
    print(f'{output[key]}')