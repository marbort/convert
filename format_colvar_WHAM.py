import plumed
import argparse
import numpy as np

    



parser = argparse.ArgumentParser(description='Format colvar file for use with WHAM')

parser.add_argument('--input',dest='input',type=str,default='COLVAR',help='full COLVAR file')
parser.add_argument('--output',dest='output',default='WHAM.dat',type=str,help='WHAM formatted file')
parser.add_argument('--hist',dest='hist',default='hist.dat',type=str,help='umbrella histogram file')
parser.add_argument('--cv', dest='cv', default="",nargs='+',
                    type=str, help='define CVs to be printed')
parser.add_argument('--split', dest='split', default=1,
                    type=int, help='number of intervals to split simulation')
parser.add_argument('--random', dest='random', 
                    action='store_true', help='randomize time series')
args = parser.parse_args()

cv=plumed.read_as_pandas(args.input)
cols=["time"]+args.cv
cv.to_csv("{}".format(args.output),index=False,columns=cols,sep=' ',header=["#time"]+args.cv, quoting=None,quotechar=' ',float_format='%10.6f')
chunk_length=len(cv)//args.split
chunks=list(range(0,len(cv)+1,chunk_length))
print(len(chunks))


if args.random:
    cv_equil=cv[chunks[1]:len(cv)]
    cv_equil.to_csv("{}".format("{}_equil".format(args.output)),index=False,columns=cols,sep=' ',header=["#time"]+args.cv, quoting=None,quotechar=' ',float_format='%10.6f')
    cv_equil_rand=cv_equil.sample(frac = 1)
    chunk_length_rand=len(cv_equil_rand)//(args.split-1)
    chunks_rand=list(range(0,len(cv_equil_rand),chunk_length_rand))
    cv[chunks[1]-chunk_length:chunks[1]].to_csv("{}_rand_{}".format(args.output,0),index=False,columns=cols,sep=' ',header=["#time"]+args.cv, quoting=None,quotechar=' ',float_format='%10.6f')
    for i,chunk in enumerate(chunks_rand[1:]):
        try:
            cv_equil_rand[chunk-chunk_length_rand:chunk].to_csv("{}_rand_{}".format(args.output,i+1),index=False,columns=cols,sep=' ',header=["#time"]+args.cv, quoting=None,quotechar=' ',float_format='%10.6f')
        except:
            print("Error!! Check File format or CV name")

else:
    for i,chunk in enumerate(chunks[1:]):
        try:
            cv[chunk-chunk_length:chunk].to_csv("{}_{}".format(args.output,i),index=False,columns=cols,sep=' ',header=["#time"]+args.cv, quoting=None,quotechar=' ',float_format='%10.6f')
        except:
            print("Error!! Check File format or CV name")

hist_all=[]
for var in args.cv:
    hist,edges=np.histogram(cv[var],bins=100,density=True)
    hist_all.append([hist,edges[:-1]])
print(len(hist_all[0][1]),len(hist_all[0][0]))
for i in range(len(args.cv)):
    with open(args.hist+'_'+args.cv[i],'w') as ofile:
        ofile.write("#{:10s} {:10s}\n".format(args.cv[i],'count'))
        for j,item in enumerate(hist_all[i][0]):
            ofile.write("{:10.6f} {:10.6f}\n".format(hist_all[i][1][j],item))
    
    
    
