# importing sys
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import argparse

# adding Folder_2/subfolder to the system path
sys.path.insert(0, '/home/marco/SHARED/GitHub/mepfinder')

from mepfinder.grid_func import GridFunc
from mepfinder.flooder import Flooder

#Use the vreco tool to generate V.final.out which contains the grid and the potential Initialize a GridFunc from V.final.out.
def extract_data(input):
    data=[]
    with open(input,'r') as ifile:
        file=np.loadtxt(ifile,unpack=True)
    
    cv1=np.unique(file[0])
    cv2=np.unique(file[1])
    val=np.array(file[2])
    free_grid=val.reshape(len(cv1),len(cv2))
    

    return(cv1,len(cv1),cv2,len(cv2),free_grid)


def plot2d(x,y,maxz,value,file,labx,laby,cmap,min_path):
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
    
    
    
    
    plt.scatter(min_path[0],min_path[1],color='white',s=10)
    plt.xlim([min(x),max(x)])
    plt.ylim([min(y),max(y)])
    plt.tight_layout()
    plt.savefig('{}_minpath.png'.format(file),format='png')

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
cv1,dim1,cv2,dim2,free_grid=extract_data(file)
surface_wrong_order=np.loadtxt(args.input)
surface_right_order=[]
for i in range(dim1):
    for j in range(dim2):
        surface_right_order.append(surface_wrong_order[i+dim2*j])
np.savetxt('b.txt',surface_right_order)
        

#print(surface[0])
gf = GridFunc.from_file('b.txt')


print(gf.shape)
print("#####")
#print(dir(gf))
print("#####")
gf.save('grid.dat')
#sf=gf.to_surface
#arr=np.array(sf)
#print(arr)
# (241, 101)
#Next, you need to specify the two points between you like to find a path.
#print(gf.points[30])
crd_p1=[]
crd_p2=[]
print(args.p1,args.p2)
for coord in args.p1:
    if coord == "None":
        crd_p1.append(None)
    else:
        crd_p1.append((float(coord)-args.tol,float(coord)+args.tol))
        
for coord in args.p2:
    if coord == "None":
        crd_p2.append(None)
    else:
        crd_p2.append((float(coord)-args.tol,float(coord)+args.tol))
        
print(crd_p1)
print(crd_p2)
p1 = gf.g_minimize(crd_p1[0],crd_p1[1])
#print(p1)
# (58, 25)

p2 = gf.g_minimize(crd_p2[0],crd_p2[1])
#print(p2)
# (87, 68)
#Initialize a Flooder based on the GridFunc

flooder = Flooder(gf)
#Find the minimum energy path connecting p1 and p2

path = flooder.flood(p1, p2)
step1=19
step2=80
#print(dir(path))
print(path.points[0])
print(path.points[-1])
#print(len(list(path)))
#print(list(path)[0])
#print(list(path)[30])
#print(list(path)[-1])
path_arr=[]
for item in list(path.points):
        path_arr.append(item)
path_arr_np=np.array(path_arr)
path_arr_np_cv1=np.array(list(zip([x[0] for x in path_arr],[x[2] for x in path_arr])))
path_arr_np_cv2=np.array(list(zip([x[1] for x in path_arr],[x[2] for x in path_arr])))
np.savetxt("path.dat",path_arr_np,fmt='%.4f')
np.savetxt("path_CV1_2_3.dat",path_arr_np_cv1[step1:step2],fmt='%.4f')
np.savetxt("path_CV2_1.dat",path_arr_np_cv2[:step1],fmt='%.4f')
np.savetxt("path_CV2_4.dat",path_arr_np_cv2[step2:],fmt='%.4f')
min_path=np.loadtxt('path.dat',unpack=True)


MAX=args.max

labx="CV1"
laby="CV2"
cmap_active='rainbow'

plot2d(cv1,cv2,MAX,free_grid,file,labx,laby,cmap_active,min_path)
    
    