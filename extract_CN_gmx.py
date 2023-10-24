import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os
import glob
from itertools import combinations
import regex as re
from scipy.signal import find_peaks
import json

def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    np.hanning, np.hamming, np.bartlett, np.blackman, np.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return(y[((window_len-1)//2):-((window_len-1)//2)])

def extract_coord_number(files):
    CN={}
    
    
    for file in files:
        split=file.split('_')[5]
        
        with open(file,'r') as ifile:
            lines=ifile.readlines()
        ref=[x.split()[-1][:-1] for x in lines if "subtitle" in x]
        sel=[x.split()[-1][1:-1] for x in lines if re.findall('s.*legend',x)]
        if CN.get(ref[0]):
            pass
        else:
            CN.update({ref[0]:{}})
        """
        if CN.get(ref[0]):
            pass
        else:
            CN[ref[0]]={}
        if CN[ref[0]].get(split):
            pass
        else:
        """
        
        
        
        x=[float(k.split()[0]) for k in lines if "#" not in k if "@" not in k ]
        rdf=[[float(k.split()[j]) for k in lines if "#" not in k if "@" not in k] 
             for j in range(1,len(sel)+1)]
        rdf_smooth=[smooth(np.array(x),11,'hanning') for x in rdf]
        mins=[find_peaks([-x for x in j],width=5) for j in rdf_smooth]
        maxs=[find_peaks(j,width=4) for j in rdf_smooth]
        print(mins[2][0][0])
        with open(file.replace('rdf','coord'),'r') as ifile:
            lines=ifile.readlines()
        coord=[[float(k.split()[j]) for k in lines if "#" not in k if "@" not in k] 
             for j in range(1,len(sel)+1)]
        #print(ref[0],sel[0],coord[0][mins[0][0][0]])
        CoordNumb=[coord[i][x[0][0]] for i,x in enumerate(mins)]
        #print(CoordNumb[0],sel[0],x[mins[0][0][0]])
        CN[ref[0]].update({split:{k:{'max':x[maxs[i][0][0]],'min':x[mins[i][0][0]],'val':CoordNumb[i]} for i,k in enumerate(sel)}})
        #CN[ref[0]][split]={x:CoordNumb[i] for i,x in enumerate(sel)}
    return(CN)

def block_avg_CN(CN):
    lst=[]
    for ref in CN:
        for sel in CN[ref][list(CN[ref].keys())[0]]:
            avg_CN=np.average([CN[ref][x][sel]['val'] for x in CN[ref]])
            std_CN=np.std([CN[ref][x][sel]['val'] for x in CN[ref]])
            lst.append((ref,sel,avg_CN,std_CN))
            
    return(lst)

def format_CN(CN):
    tot_split=len(CN[list(CN.keys())[0]])
    with open('CN_summary.csv','w') as ofile:
        ofile.write("Ref,Sel,")
        for i in range(tot_split):
            ofile.write("Split{}_min,".format(i))
        for i in range(tot_split):
            ofile.write("Split{}_max,".format(i))
        for i in range(tot_split):
            ofile.write("Split{}_CN,".format(i))
        ofile.write("\n")
        for ref in CN:
            for sel in CN[ref][list(CN[ref].keys())[0]]:
                ofile.write("{},{},".format(ref,sel))
                for split in CN[ref]:
                    ofile.write("{},".format(CN[ref][split][sel]['min']))
                for split in CN[ref]:
                    ofile.write("{},".format(CN[ref][split][sel]['max']))
                for split in CN[ref]:
                     ofile.write("{},".format(CN[ref][split][sel]['val']))
                ofile.write("\n")
                    
                
            
        
        





CN=extract_coord_number(glob.glob('rdf*all*.xvg'))
#print(CN)
lst=block_avg_CN(CN)
with open('coords_split.json','w') as ofile:
    json.dump(CN,ofile,indent=2)
format_CN(CN)
print(lst)

