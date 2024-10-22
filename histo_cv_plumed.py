import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import plumed
import argparse
import sys
import numpy as np


data=np.loadtxt(sys.argv[1],unpack=True)

histo,bins=np.histogram(data[int(sys.argv[21])],bins=4,range=(0,4))

fig=plt.figure(figsize=(16,12),dpi=150)
font = {'family' : 'Formular',
    'weight' : 'normal',
    'size'   : 46}
mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['lines.linewidth'] = 3

plt.bar(bins[1:],histo/len(data[9]),width=0.5)
plt.tight_layout()
plt.show()
