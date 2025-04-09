import sys
import numpy as np
import shutil

shutil.copy(sys.argv[1],sys.argv[1]+'.bak')

with open(sys.argv[1], 'r') as f:
    lines = f.readlines()

col=lines[0].split().index("CV2")-2
print(col)
newlines    = []
for line in lines:
    if line.startswith('#'):
        newlines.append(line)
    elif float(line.split()[col]) < 1.5:
        pass
    else:
        newlines.append(line)
with open(sys.argv[1], 'w') as f:
    for line in newlines:
        f.write(line)

    
