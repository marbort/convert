import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import plumed
import argparse


def extract_data(input):
    colvar=plumed.read_as_pandas(input)
    return(colvar)

def calculate_2d_hist(data,cvs,temp):
    beta=1./(0.00831441001626*temp)
    logweights=beta*data['opes.bias']
    logweights -= np.amax(logweights)
    histo, bin_edges_cv1, bin_edges_cv2 = np.histogram2d(data[cvs[0]],data[cvs[1]],weights=np.exp(logweights),bins=50)
    #err = np.sqrt(np.histogram2d(data[cvs[0]],data[cvs[1]],weights=np.power(np.exp(logweights),2),bins=50))
    bin_centers_cv1 = (bin_edges_cv1[1:]+bin_edges_cv1[:-1])/2
    bin_centers_cv2 = (bin_edges_cv2[1:]+bin_edges_cv2[:-1])/2
    fes = -(1/beta)*np.log(histo)
    offset = np.min(np.ma.masked_invalid(fes))
    fes -= offset
    feserr = 0 #(1/beta)*err/histo
    np.savetxt('fes.dat',fes,fmt='%.5f')
    np.savetxt('bins_{}.dat'.format(cvs[0]),bin_centers_cv1,fmt='%.2f')
    np.savetxt('bins_{}.dat'.format(cvs[1]),bin_centers_cv2,fmt='%.2f')
    return(fes,feserr,bin_centers_cv1,bin_centers_cv2)


def plot2d(x,y,value):
    fig=plt.figure(figsize=(16,10),dpi=150)
    lev=int(round(np.max(np.ma.masked_invalid(value))/4,0))
    #plt.imshow(value,extent=(min(x),max(x),min(y),max(y)))
    plt.contourf(x, y,value,lev)
    plt.colorbar()
    plt.savefig('fes.png',format='png')



def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--input' , dest='input',help='string to match input files')
    parser.add_argument('--plot' , dest='plot', action='store_true',help='plot ')
    parser.add_argument('--cvs' , dest='cvs', nargs='+' , help='cvs to plot')
    parser.add_argument('--temp' , dest='temp', help='temp (K)',default=300,type=float)
    parser.add_argument('--wd' , dest='wd', help='plot width',default=None)
    parser.add_argument('--hg' , dest='hg', help='plot height',default=None)
    parser.add_argument('--fac' , dest='fac', help='CV conversion factor',default=1, type=float)
    args = parser.parse_args()
    
    data=extract_data(args.input)
    fes,feserr,bin_centers_cv1,bin_centers_cv2=calculate_2d_hist(data,args.cvs,args.temp)
    plot2d(bin_centers_cv1,bin_centers_cv2,fes)
    
    #print(bin_centers)
    

if __name__=="__main__":
    main()
