import os
import numpy as np
import sys


with open(sys.argv[1],'r') as ifile:
    lines=ifile.readlines()
if sys.argv[2]=="angle":
    x=[float(x.split()[0])*0.0174533 for x in lines]
else:
    x=[float(x.split()[0]) for x in lines]
y=[float(x.split()[1]) for x in lines]
    
fit=np.polyfit(x,y,2)
print(fit)
    
    
