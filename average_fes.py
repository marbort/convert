import numpy as np
import argparse

def average_fes(paths,oned=False):
    fes = []
    std_fes = []
    for path in paths:
        print(path)
        fes.append(np.loadtxt(path))
    fes = np.array(fes)
    avg_fes = np.mean(fes, axis=0)
    if oned:
        std_fes = np.std(fes, axis=0)
    return avg_fes,std_fes
#symmetrize the FES respect to the x=y line
def symmetrize_fes(fes):
    cv1_temp = fes[:,0]
    cv2_temp = fes[:,1]
    cv1 = np.unique(cv1_temp)
    cv2 = np.unique(cv2_temp)
    val = fes[:,2]
    free_grid = val.reshape(len(cv1), len(cv2))
    print(free_grid.shape)
    free_grid_sym = np.zeros((len(cv1), len(cv2)))
    for i in range(len(cv2)):
        for j in range(len(cv1)):
            free_grid_sym[i,j] = 0.5*(free_grid[i,j] + free_grid[j,i])
    fes_sym = np.vstack((cv1_temp,cv2_temp,free_grid_sym.flatten()))
    return fes_sym.T

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('fes', nargs='+', help='FES files to average')
    parser.add_argument('-o', '--output', help='Output file')
    parser.add_argument('--symm', action='store_true', help='Symmetrize the FES')
    parser.add_argument('--oned', action='store_true', help='1D FES')
    
    args = parser.parse_args()

    avg_fes,std_fes = average_fes(args.fes,args.oned)
    if args.oned:
        avg_fes[:,-1] = std_fes[:,-2]
    
    if args.symm:
        avg_fes_sym = symmetrize_fes(avg_fes)
    if args.output:
        np.savetxt(args.output, avg_fes,fmt='%10.6f')
        if args.symm:
            np.savetxt(args.output + '_sym', avg_fes_sym,fmt='%10.6f')
    else:
        print(avg_fes)

if __name__ == '__main__':
    main()  