import numpy as np
import argparse
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
from matplotlib import font_manager
import json
import sys

with open(sys.argv[1], 'r') as f:
    data = json.load(f)

err_MC = np.loadtxt(sys.argv[2],unpack=True)[2]
print(err_MC)

print(np.subtract(data['summary']['equil']['free'],data['summary']['std_free_mean'])[0])

fig=plt.figure(figsize=(10,10))
plt.plot(data['summary']['equil']['cv'],data['summary']['equil']['free'],label='Free Energy',color='b')
plt.fill_between(data['summary']['equil']['cv'], np.add(data['summary']['equil']['free'],data['summary']['std_free_mean'])                        ,
                np.subtract(data['summary']['equil']['free'],data['summary']['std_free_mean']),
                linewidth=0,
                alpha=0.1,
                color='b')
plt.fill_between(data['summary']['equil']['cv'], np.add(data['summary']['equil']['free'],err_MC)                        ,
                np.subtract(data['summary']['equil']['free'],err_MC),
                linewidth=0,
                alpha=0.1,
                color='r')
                    
plt.show()

