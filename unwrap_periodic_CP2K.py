import numpy as np
import os
import argparse

def unwrap_periodic(paths, threshold,outdir,pos,neg):
    for file in paths:
        data=np.loadtxt(file,unpack=True)
        for i,val in enumerate(data[1]):
            if neg:
                if val > threshold:
                    data[1][i] -= 2*np.pi
            else:
                if val < threshold:
                    data[1][i] += 2*np.pi
        if pos is not None:
            outdir=os.path.join(outdir,str(data[0][pos]))
            val=int(os.path.basename(file).split('_')[pos])
            name=os.path.basename(file).replace(val,f'{val:3.d,}')
            np.savetxt(os.path.join(outdir,name),data.T,fmt='%15.8f')
        np.savetxt(os.path.join(outdir,os.path.basename(file)),data.T,fmt='%15.8f')


def main():
    
    parser = argparse.ArgumentParser(description='Unwrap periodic data')
    parser.add_argument('paths', type=str, nargs='+', help='Path to the files to unwrap')
    parser.add_argument('-t', '--threshold', type=float, default=0, help='Threshold for unwrapping')
    parser.add_argument('--neg', action='store_true', help='Unwrap positive values')
    parser.add_argument('-o', '--outdir', type=str, default='.', help='Output directory')
    parser.add_argument('--pos',type=int, default=None, help='Position of the angle value in the file name.')
    
    args=parser.parse_args()
    
    unwrap_periodic(args.paths,args.threshold,args.outdir,args.pos,args.neg)
    



if __name__ == "__main__":
    main()
            
            
    