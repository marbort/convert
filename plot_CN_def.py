import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import argparse
import sys
import numpy as np


font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 22}
mpl.rc('font', **font)

def coord_number(d,d0,p,q):
    CN=(1-(d/d0)**p)/(1-(d/d0)**q)
    return(CN)

def plot_CN(dists,val,label_name,label_val):
    plt.plot(dists,val,label="{}={}".format(label_name,label_val))
    
    
p=[6,12,24,48]
q=[x*2 for x in p]

d0=range(1,5)


fig=plt.figure(figsize=(16,10),dpi=150)
for i,item in enumerate(p):
    dists=np.linspace(0,d0[0]+1,200)
    CN=[coord_number(x,d0[0],item,q[i]) for x in dists]
    plot_CN(dists,CN,"p",item)
plt.xlabel("Bond Distance ($\AA$)")
plt.ylabel("Coordination Number")
plt.legend()
plt.title("$d_0=1.0\ \AA$")
plt.savefig('CN_p.png',format='png')
fig=plt.figure(figsize=(16,10),dpi=150)
for i in d0:
    dists=np.linspace(0,i+2,200)
    CN=[coord_number(x,i,p[1],q[1]) for x in dists]
    plot_CN(dists,CN,"$d_0$",i)
plt.title("p=12, q=24")
plt.xlabel("Bond Distance ($\AA$)")
plt.ylabel("Coordination Number")
plt.legend()
plt.savefig('CN_d0.png',format='png')
    
    



    
