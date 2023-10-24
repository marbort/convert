#%%
import numpy as np
import matplotlib.pyplot as plt
import os


root="/home/marco/SHARED/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/IMC/umbrella_30_200/"
const=200
start=9.13611
with open(root+"window2/umbrella_pullf.xvg") as ifile:
    lines=ifile.readlines()
t=[float(x.split()[0]) for x in lines if "@" not in x if "#" not in x]
f=[float(x.split()[1]) for x in lines if "@" not in x if "#" not in x]

x=[k/const+start for k in f]

hist=np.histogram(x,200,(8.48,10.00))

plt.plot(hist[1][:-1],hist[0])

with open(root+"hist_test.xvg") as ifile:
    lines=ifile.readlines()
    hist_x=[float(x.split()[0]) for x in lines if "@" not in x if "#" not in x]
    hist_y=[float(x.split()[2]) for x in lines if "@" not in x if "#" not in x]

plt.plot(hist_x,hist_y,label="WHAM")
plt.legend()

# %%
