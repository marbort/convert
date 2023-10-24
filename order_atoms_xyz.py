import numpy as np
import sys

def reorder(input,ref):
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    data=lines[2:]
    #print(data)
    ref=np.array([float(data[ref].split()[x]) for x in range(1,4)])
    dist=[(np.linalg.norm(np.array([(float(x.split()[k])) for k in range(1,4)]) -ref),i) for i,x in enumerate(data)]
    dist.sort()
    return(data,dist)
def print_ordered(data,order,output):
    with open(output,'w') as ofile:
        ofile.write("{}\n".format(len(data)))
        ofile.write("Ordered Atoms\n")
        for line in order:
            ofile.write("{}".format("".join(data[line[1]])))


def main():
    data,order=reorder(sys.argv[1],81)
    print_ordered(data,order,sys.argv[2])


if __name__=="__main__":
    main()

    