import math

with open('colvar','r') as ifile:
    lines=ifile.readlines()
print(math.floor((len(lines)-1000)/10))