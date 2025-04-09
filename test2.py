import sys

with open('paths.txt','r') as ifile:
    paths = ifile.readlines()
    print(paths)
    for line in paths:
        print(line.strip())
        