import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

def extract_data(input):
    data={}
    filename=os.path.splitext(input)[0]
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    for i in range(len(lines[0].split())):
            data[i]=[float(x.split()[i]) for x in lines]
    print(data)
    return(data)


def plot_histo(data):
    hist_all=[]
    for i in range(1,len(data)):
        hist,edges=np.histogram(data[i],bins=50)
        hist_all.append([hist,edges[:-1]])
    print(hist_all)
    fig=plt.figure()
    for i in hist_all:
        plt.plot(i[1],i[0])
    plt.show()
        
        

def main():
    parser = argparse.ArgumentParser(description='Plot data')


    parser.add_argument('--input', dest='input', default=None,
                    type=str, help='input file')

    parser.add_argument('--xcols', dest='xcols', default=None,nargs='+',
                    type=int, help='x columns')
    parser.add_argument('--ycols', dest='ycols', default=None,nargs='+',
                    type=int, help='y columns')
    args = parser.parse_args()
    
    data=extract_data(args.input)
    plot_histo(data)

if __name__ == "__main__":
   main()