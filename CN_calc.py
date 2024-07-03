import numpy as np
import sys


def coord_Number(x,y,dist_min,dist_max):
    val_x=[k for k in x if dist_min<= k <=dist_max]
    val_x_sq=[k**2 for k in val_x]
    val_y=[y[i] for i,k  in enumerate(x) if dist_min <= k <= dist_max]
    prod=[k*val_y[i] for i,k in enumerate(val_x_sq)]
    CN_all=4*np.pi*np.trapz(prod,val_x)
    return(CN_all)


def extract_data_xmgrace(input):
    with open(input,'r') as ifile:
        lines=ifile.readlines()

    # Generate example data (you can replace this with your own data)
    x = [float(x.split()[0]) for x in lines if ("&" not in x) and ("@" not in x) and ("#" not in x)]
    y = [float(x.split()[1]) for x in lines if ("&" not in x) and  ("@" not in x) and ("#" not in x)]
    
    return(x,y)

def main():
    input=sys.argv[1]
    min_loc=float(sys.argv[2])

    
    x,y=extract_data_xmgrace(input)
    CN=coord_Number(x,y,0,min_loc)
    print(CN)
    
if __name__=="__main__":
    main()

    




