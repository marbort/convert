import argparse
import numpy as np

    



parser = argparse.ArgumentParser(description='Format colvar file for use with WHAM')

parser.add_argument('--input',dest='input',type=str,default='pullx.xvg',help='pullx file')
parser.add_argument('--output',dest='output',default='WHAM.dat',type=str,help='WHAM formatted file')
parser.add_argument('--hist',dest='hist',default='hist.dat',type=str,help='umbrella histogram file')
parser.add_argument('--cv', dest='cv', default="",nargs='+',
                    type=str, help='define CVs to be printed')
parser.add_argument('--split', dest='split', default=1,
                    type=int, help='number of intervals to split simulation')
parser.add_argument('--random', dest='random', 
                    action='store_true', help='randomize time series')
args = parser.parse_args()

cv=np.loadtxt(args.input,comments=["#","@"])
chunk_length=len(cv)//args.split
chunks=list(range(0,len(cv)+1,chunk_length))
print(len(chunks))
np.savetxt("{}".format(args.output),cv,fmt='%.6f')


for i,chunk in enumerate(chunks[1:]):
        try:
            np.savetxt("{}_{}".format(args.output,i),cv[chunk-chunk_length:chunk],fmt='%.6f')
        except:
            print("Error!! Check File format or CV name")


    
