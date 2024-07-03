import pandas as pd
import numpy as np
import glob
import sys
import matplotlib.pyplot as plt

def rational_function(range:tuple,r0,NN,MM):
    x=np.linspace(min(range),max(range),100)
    y=np.fromiter((((1-(xi/r0)**NN)/(1-(xi/r0)**MM)) for xi in x),x.dtype)
    return(x,y)


def main():
    r0=np.linspace(2.5,3,6)
    NN=[6,12,24]
    MM=[2*x for x in NN]
    print(MM)

    inputs=glob.glob(sys.argv[1])
    maxim=0
    minim=100
    fig=plt.figure(figsize=(16,10),dpi=150)
    
    for file in inputs:
        df=pd.read_csv(file,sep=" ",header=None)
        print(df[0])
        if df.max()[0] > maxim:
            maxim=df.max()[0]
        if df.min()[0] < minim:
            minim=df.min()[0]
        plt.plot(df[0],df[1]/df[1].max())
    print(minim,maxim)
    limits=(minim-1,maxim+1)
    rationals=[(rational_function(limits,x,y,MM[j]),x,y,MM[j]) for x in r0 for j,y in enumerate(NN)]
    
    for rat in rationals:
        plt.plot(rat[0][0],rat[0][1],label=f"r0={rat[1]},NN={rat[2]},MM={rat[3]}")
        plt.legend()
    
    plt.show()
    

if __name__ == "__main__":
    main()






    


