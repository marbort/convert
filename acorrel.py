import numpy as np
import argparse 


def extract_data(input,col):
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    data=[float(x.split()[col]) for x in lines if '#' not in x]
    return(data)

def acorrel(data):

    x = np.array(data) 

    # Mean
    mean = np.mean(data)

    # Variance
    var = np.var(data)

    # Normalized data
    ndata = data - mean

    acorr = np.correlate(ndata, ndata, 'full')[len(ndata)-1:] 
    acorr = acorr / var / len(ndata)
    return(acorr)

def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--input', dest='input', 
                        type=str, help='input data')
    parser.add_argument('--col', dest='col', default=1,
                        type=int, help='column storing time series')
    args = parser.parse_args()
    
    data=extract_data(args.input,args.col)
    acorr=acorrel(data)
    with open('autocoor.dat','w') as ofile:
        for line in acorr:
            ofile.write("{}\n".format(line))
    #print(acorr)
    
if __name__=="__main__":
    main()
    

    
    
 