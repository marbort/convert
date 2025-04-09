import numpy as np
import argparse



def scale_convert_plot_wham(file,min,factorx=1,factory=1):
    data=np.loadtxt(file,comments="#",unpack=True)
    data[0]=data[0]*factorx
    if min is not None:
        diffs=abs(data[0]-min)
        min_diffs_idx=np.argmin(diffs)
        min_value=data[1][min_diffs_idx]
        data[1]=data[1]-min_value
    else:
        data[1]=data[1]-np.min
    data[1]=data[1]*factory
    return data
    

def calculated_error(file):
    data=[]
    for i in range(1,100):
        try:
            data.append(np.loadtxt(f"{file}_{i}",comments="#",unpack=True)[1])
        except:
            print(f"Found {i-1} files")
            break
    error=(np.std(data,axis=0)/np.sqrt(len(data)))
    return error
        
    

def main():
    parser = argparse.ArgumentParser(description='Scale and convert WHAM data')
    parser.add_argument('data', type=str, help='Path to the data file')
    parser.add_argument('--min', type=float, help='Minimum value',default=None)
    parser.add_argument('--factorx', type=float, default=1, help='Factor to scale X data')
    parser.add_argument('--factory', type=float, default=1, help='Factor to scale Y data')
    parser.add_argument('--splits', type=str, help='Splits filename', default=None)
    
    args = parser.parse_args()
    
    data=scale_convert_plot_wham(args.data,args.min,args.factorx,args.factory)
    error = calculated_error(args.splits)
    data = np.vstack((data,error))
    np.savetxt(f"scaled_{args.data}", data.T, fmt='%10.6f')
    
if __name__ == "__main__":
    main()
    
    
    
    
    