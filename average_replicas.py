import numpy as np
import glob
import sys


inputs=glob.glob(sys.argv[1])
print(inputs)
data=[]
for input in inputs:
    data.append(np.loadtxt(input,unpack=True))
#print(data[4])
avg=np.average(data,axis=0)
print(avg[-1])
avg[-1]=avg[-1]-min(avg[-1])
np.savetxt(sys.argv[2],np.transpose(avg),fmt='%.3f')


        
        
    
    


