import numpy as np
import sys

def find_first_max_index(array):
    value = -1000
    for val in array:
        #print(val)
        if val >= value:
            #print(val)
            value = val
        else:
            break
    return (value,array.tolist().index(value))

def find_first_min_index(array, maxindex):
    value = 1000
    for val in array[maxindex:]:
        if float(val) <= value:
            value = val
        else:
            break
    return (value,array.tolist().index(value))


file=np.loadtxt(sys.argv[1], comments=['#','@'],unpack=True)
#print(file[1])
max,maxindex = find_first_max_index(file[1])
min,minindex = find_first_min_index(file[1],maxindex)

file2 = np.loadtxt(sys.argv[2], comments=['#','@'],unpack=True)

#print(f"File,Max,Min,CN")
print(f"{sys.argv[1]},{sys.argv[2]},{file[0][maxindex]},{file[0][minindex]},{file2[1][minindex-1]}")




