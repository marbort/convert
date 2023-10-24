import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import regex as re

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
    return( y[((window_len-1)//2):-((window_len-1)//2)])


file=sys.argv[1]
with open(file,'r') as ifile:
            out_name=ifile.name.replace('.xvg','_smooth.xvg')
            lines=ifile.readlines()
sel=[x.split()[-1][1:-1] for x in lines if re.findall('s.*legend',x)]
parms=[k for k in lines if "#" in k or "@" in k  ]
x=[float(k.split()[0]) for k in lines if "#" not in k if "@" not in k ]
rdfs=[[float(k.split()[j]) for k in lines if "#" not in k if "@" not in k] 
             for j in range(1,len(sel)+1)]
windows=['flat', 'hanning', 'hamming', 'bartlett', 'blackman']
rdf_conv=[smooth(np.array(rdf),int(sys.argv[2]),'hanning') for rdf in rdfs]
xn=np.linspace(min(x),max(x),len(rdf_conv[0]))
#print(len(rdf))
with open(out_name,'w') as ofile:
    for parm in parms:
        ofile.write(parm)
    for i,crd in enumerate(x):
        ofile.write("{}".format(crd))
        for k in range(len(sel)):
            ofile.write("    {}".format(rdf_conv[k][i]))
        ofile.write("\n")
"""
for i,graph in enumerate(rdf_conv):
    print(len(rdf_conv))
    plt.plot(x,graph,label=windows[i])
plt.plot(x,rdf,label='orig')
plt.legend()
plt.show()
"""
