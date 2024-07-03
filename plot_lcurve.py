
import argparse
import numpy as np
import matplotlib.pyplot as plt
import glob
import pandas as pd


data=np.loadtxt('lcurve.out',unpack=True)


print(data)


fig=plt.figure(figsize=(16,9),dpi=150)
labels=["rmse_trn","rmse_e_trn","rmse_f_trn","lr"]
limits=[[-1,10],[-0.1,0.5],[-0.1,0.5],[-1e-5,1e-3]]
for i in range(1,5):
    plt.subplot(2,2,i)
    plt.plot(data[0],data[i])
    plt.title(labels[i-1])
    plt.ylim(limits[i-1])
plt.show()




    





        

    









    
    
    
