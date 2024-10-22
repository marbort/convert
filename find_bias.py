import sys

def find_bias(input):
    with open(input,'r') as ifile:
        header=ifile.readline()
    bias=",".join([x for x in header.split() if "bias" in x])
    return(bias)

bias=find_bias('colvar')
print(bias)