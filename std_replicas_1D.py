import numpy as np
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl




inputs=glob.glob(sys.argv[1])
print(inputs)
data=[]
for input in inputs:
    data.append(np.loadtxt(input,unpack=True))
#print(data[4])
avg=np.average(data,axis=0)
dev=np.std(data,axis=0)

np.savetxt(f'{sys.argv[2]}_st_dev.dat',np.transpose([data[0][0],dev[1]]),fmt='%.3f')


fig=plt.figure(figsize=(16,10),dpi=150)
font = {'family' : 'sans',
    'weight' : 'normal',
    'size'   : 32}
mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 3

plt.plot(avg[0],avg[1],label='Average')
plt.fill_between(avg[0],avg[1]-dev[1],avg[1]+dev[1],alpha=0.3)
for i,item in enumerate(data):
    plt.plot(item[0],item[1],alpha=0.5,linewidth=0.5,label=str(i))
plt.ylim([-10,200])

plt.legend()
plt.savefig(f"{sys.argv[2]}_std.png",format='png',dpi=150)

        
        
    
    


