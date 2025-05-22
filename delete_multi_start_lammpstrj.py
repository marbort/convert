import sys
import regex as re

begin=0

with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
match = re.compile(r'^0\n')
for i,line in enumerate(lines):
        match.search(line)
        if match.search(line):
            print("Found")
            begin=i-1

with open(sys.argv[1], 'w') as of:
    for line in lines[begin:]:
        of.write(line)
