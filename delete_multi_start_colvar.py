import sys
import regex as re

begin=1

with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
len_line = len(lines[2].split())
print(len(lines[0].split()),len(lines[1].split()))
match = re.compile(r'^\ 0.000000')
for i,line in enumerate(lines):
    if i == 0:
        continue
    if len(line.split()) != len_line:
            del lines[i]
    else:
        match.search(line)
        if match.search(line):
            begin=i-1

with open(sys.argv[1], 'w') as of:
    for line in lines[begin:]:
        of.write(line)
