from datetime import datetime
import argparse
from timediff import time_diff


def extract_exc_time(inputs):
    starts=[]
    ends=[]
    files=[]
    for input in inputs:
        with open(input,'r') as ifile:
            lines=ifile.readlines()
        for line in lines:
            if "PROGRAM STARTED AT" in line:
                start=line.split()[-2]+" "+line.split()[-1][:-4]
            if "PROGRAM ENDED AT" in line:
                end=line.split()[-2]+" "+line.split()[-1][:-4]
        files.append(input)
        starts.append(start)
        ends.append(end)
    return(files,starts,ends)
            



def main():
    parser = argparse.ArgumentParser(
                    prog='CP2K_time',
                    description='Calculate execution time in seconds',
                    epilog='Time is nothing, timing is everything')
    parser.add_argument('--input', help='CP2K output files', nargs='+')
    args = parser.parse_args()
    files,starts,ends=extract_exc_time(args.input)
    for i,item in enumerate(starts):
        t1,t2,delta=time_diff(item,ends[i])
        print('CP2K output: ', files[i])
        print('Start time:', t1)
        print('End time:', t2)
        print(f"Time difference is {delta.total_seconds()} seconds\n")
    

    

if __name__=="__main__":
    main()

