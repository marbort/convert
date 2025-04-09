import math
import sys
import numpy as np
import os


def get_data(paths,cvmax):
    data={}
    inputs=[]
    with open(paths,'r') as ifile:
        lines=ifile.readlines()
        for line in lines:
            inputs.append(line.strip())
    for input in inputs:
        with open(input,'r') as ifile2:
            lines=ifile2.readlines()
        for j,line in enumerate(lines):
            if "Composition" in line:
                try:
                    for i in range(int(cvmax)+1):
                        data[line.split('at')[-1].strip()].append(lines[j+i+1].split(':')[-1])
                    print(f"Found {line.split('at')[-1].strip()} in {input}")
                except:
                    try:
                    
                        data[line.split('at')[-1].strip()]=[]
                        for i in range(int(cvmax)+1):
                            data[line.split('at')[-1].strip()].append(lines[j+i+1].split(':')[-1])
                    except:
                        print(f"Error in {input} with {line.split('at')[-1].strip()}")
    #print(data)
    return data

def get_mean_data(data,cvmax,files):
    mean_data={}
    std_data={}
    #print(data)
    for key in data.keys():
        mean_data[key]=np.mean(np.array(data[key],dtype=float).reshape(files,cvmax+1),axis=0)
        std_data[key]=np.std(np.array(data[key],dtype=float).reshape(files,cvmax+1),axis=0)
    #print(data)
    #print(mean_data)
    return mean_data,std_data

def write_data(output,mean_data,std_data,cvmax,basins):
        for key in mean_data.keys():
            output.write(f"Composition of minimum at {key}\n")
            for j in range(cvmax+1):
                output.write(f"CN {j:.2f}: {mean_data[key][j]:.3f} {std_data[key][j]:.3f} \n")
            output.write("\n")

            
def main():
    paths=sys.argv[1]
    cvmax=sys.argv[2]
    basins=int(sys.argv[3])
    with open (paths,'r') as ifile:
        lines=ifile.readlines()
        files=len(lines)
    #files=int(sys.argv[4])
    data=get_data(paths,cvmax)
    mean_data,std_data=get_mean_data(data,int(cvmax),files)
    output=open("basin_analysis.txt",'w')
    write_data(output,mean_data,std_data,int(cvmax),basins)
    output.close()

if __name__ == "__main__":
    main()



    
    
  