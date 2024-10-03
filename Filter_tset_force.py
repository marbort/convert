import numpy as np
import dpdata as dp
import glob
import sys
import argparse



def filter_tset(input,max_force):
    unfiltered_dpdata=dp.LabeledSystem(input,'deepmd/npy')
    filtered_indices = []
    for i in range(len(unfiltered_dpdata['forces'])):
        if np.linalg.norm(unfiltered_dpdata['forces'][i], axis=1).max() < max_force:
            filtered_indices.append(i)
    print(f"Total number of configurations in tset: {len(unfiltered_dpdata)}")
    if len(filtered_indices) == 0:
        print(f"No configurations with forces less than {max_force} found in {input}")
    else:
        print(f"{len(filtered_indices)} configurations with forces less than {max_force} found in {input}")
        filtered_dpdata = unfiltered_dpdata.sub_system(filtered_indices)
        filtered_dpdata.to('deepmd/npy', f'Filtered_{input}')
        

    
    


def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--input', dest='input', 
                        type=str, help='input data',nargs="+")
    parser.add_argument('--max_force', dest='max_force', 
                        type=float, help='max force')
    args = parser.parse_args()
    
    
    max_force=args.max_force
    
    for input in args.input:
        files=glob.glob(input)
        for file in files:
            filter_tset(file,max_force)


if __name__ == "__main__":
    main()
        
