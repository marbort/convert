from timediff import time_diff
import regex as re
import sys

def read_output(output):
    with open(output,'r') as ifile:
        lines=ifile.read()
    start= " ".join([re.search('STARTED AT.*',lines)[0].split()[-2],re.search('STARTED AT.*',lines)[0].split()[-1][:-4]])
    end  = " ".join([re.search('ENDED AT.*',lines)[0].split()[-2],re.search('ENDED AT.*',lines)[0].split()[-1][:-4]])
    return(start,end)

def main():
    start,end=read_output(sys.argv[1])
    time=time_diff(start,end)
    
    print(f"""
        Started at {time[0]} 
        Ended   at {time[1]} 
        Total time {time[2].total_seconds()} s
        """
    )

if __name__=="__main__":
    main()
    
    
    

    