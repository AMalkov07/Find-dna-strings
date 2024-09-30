dict = {"s1": "cat",
     "s2": "dog",
     "s3": "cow"}


class customType:
    def __init__(self, input_s_li, input_s_ri, dict_str, is_reverse):
        self.input_s_li = input_s_li
        self.input_s_ri = input_s_ri
        self.dict_str = dict_str
        self.is_reverse = is_reverse
    
    def __repr__(self):
        return f"CustomType(input_s_li={self.input_s_li}, input_s_ri={self.input_s_ri}, dict_str={self.dict_str}, is_reverse={self.is_reverse})"
    
#var1 = customType(420, 69)

#print(var1)

input_s = "catdog"

max_length = 3
output = {}

def rec(input_s, offset):
    print(f"running rec on {input_s}")
    for i in range(max_length, -1, -1):
        sub_str =input_s[0:i]
        for key in dict.keys():
            if dict[key]== sub_str:
                output[sub_str] = customType(offset, offset+i, key, False)
                rec(input_s[i:], i)
                return


rec(input_s, 0)
for key in output.keys():
    print(f'{key}: {output[key]}')