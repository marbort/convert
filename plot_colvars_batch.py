import glob
import matplotlib.pyplot as plt
import sys
import numpy as np

max_it=sys.argv[1]
space=np.linspace(0,1,int(max_it))



def plot_colvar(file):
    iteration=int(file.split('_')[2])-1
    factor=
    try:
        data=np.loadtxt(file,unpack=True)
        plt.scatter(data[0],data[1],s=5,color=plt.cm.rainbow(space[iteration]),label=iteration)
    except:
        print(f"Error in file: {file}")


files=glob.glob("colvar*")
fig=plt.figure(figsize=(16,10),dpi=150)
for file in files:
    plot_colvar(file)
plt.xlim([0,500])
plt.legend(ncols=5)
plt.show()