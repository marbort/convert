import numpy as np
import sys

def reorder(input,ref):
    data=np.loadtxt(input)
    ref=[data[ref][x] for x in range(3)]
    dist=[(np.linalg.norm[x-ref],i) for i,x in enumerate(data)]
    dist.sort()
    return(data,dist)
def print_ordered(data,order,output):
    with open(output,'w') as ofile:
        ofile.write("{}\n".format(len(data)))
        ofile.write("Ordered Atoms\n")
        for line in order:
            ofile.write("{}\n".format("  ".join(data[line[1]])))


def main():
    data,order=reorder(sys.argv[1])
    print_ordered(data,order,sys.argv[2])


if __name__=="__main__":
    main()

    