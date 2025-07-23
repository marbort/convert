import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

def extract_data(input):
    data=[]
    with open(input,'r') as ifile:
        file=np.loadtxt(ifile,unpack=True)
    
    cv1=np.unique(file[0])
    cv2=np.unique(file[1])
    val=np.array(file[2])
    free_grid=val.reshape(len(cv1),len(cv2))
    

    return(cv1,len(cv1),cv2,len(cv2),free_grid)

def plot2d(x,y,maxz,value,file,labx,laby,cmap):
    fig=plt.figure(figsize=(16,12),dpi=150)
    font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 46}
    mpl.rc('font', **font)
    mpl.rcParams['axes.linewidth'] = 3
    mpl.rcParams['lines.linewidth'] = 3
    #lev=int(round(np.max(np.ma.masked_invalid(value))/10,0))
    MAX=int(maxz)
    
    
    lev=range(0,MAX+5,5)
    
    value[value>MAX]=MAX+10
    CLines=plt.contour(x, y,value,levels=range(0,MAX,20),vmin=0,vmax=MAX,linewidths=1.5,colors='black')
    
    print(max(x),max(y))
    tix=np.linspace(0,max(x),round(max(x)/0.5)+1)
    tiy=np.linspace(0,max(y),round(max(y)/0.5)+1)
    print(tix)
   
    plt.contourf(x, y,value,lev,vmin=0,vmax=MAX,cmap=cmap)
  
    
    plt.xlabel(labx)
    plt.ylabel(laby)
    
    
    cbar=plt.colorbar(label="$\Delta A\ (kJ\ mol^{-1})$",ticks=range(0,MAX+20,20))
    #cbar.ax.set_ylim(0,MAX)
    
    plt.xlim([min(x),max(x)])
    plt.ylim([min(y),max(y)])
    plt.tight_layout()
    plt.savefig(f"{file.replace('.dat','.png')}",format='png')

def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--input', dest='input', 
                        type=str, help='input data')
    parser.add_argument('--output', dest='output', 
                        type=str, help='output file',default="2D_reduced_compare.png")
    parser.add_argument('--title', dest='title', 
                        type=str, help='plot title',default="")
    parser.add_argument('--labels', dest='labels', default=None, type=str, nargs='+',
                            help='Plot labels label')
    parser.add_argument('--xlab', dest='xlab', type=str, 
                            help='X axis label',default='CV1')
    parser.add_argument('--ylab', dest='ylab', type=str, 
                            help='Y axis label',default='CV2')
    parser.add_argument('--min', dest='min', 
                        type=float, help='CV value for 0')
    parser.add_argument('--limy', dest='limy', 
                        type=float, help='y MAX value')
    parser.add_argument('--max', dest='max', 
                        type=int, help='z MAX value',default=200)
    parser.add_argument('--errors', dest='errors', 
                        help='plot errors',action='store_true',default=False)
    parser.add_argument('--p1', dest='p1', 
                        help='p1 coords',nargs='+',type=str)
    parser.add_argument('--p2', dest='p2', 
                        help='p2 coords',nargs='+',type=str)
    parser.add_argument('--tol', dest='tol', 
                        help='tolearnce to find min around point',type=float,default=0.1)
        
    args = parser.parse_args()
    file=args.input
    cmap_active='rainbow'
    labx=args.xlab
    laby=args.ylab
    MAX=args.max
    cv1,dim1,cv2,dim2,free_grid=extract_data(file)
    plot2d(cv1,cv2,MAX,free_grid,file,labx,laby,cmap_active)
    
if __name__ == "__main__":
    main()