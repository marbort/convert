#%%
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys

#root='/home/marco/SHARED/RATIO/WP1/ML/fox/iPrMgCl/SAGA/META/opes/Mg1-O_Mg2-O/30_it24'
#slice=sys.argv[1]
cv=2
size={}
try:
    file='fes-rew_square_sparse.dat'
    with open(file,'r') as ifile:
        lines=ifile.readlines()
        for i in range(10):
            if ' SET nbins_' in lines[i]:
                size[lines[i].split()[-2].split('_')[1]]=int(lines[i].split()[-1].rstrip())
        print(size)
except:
    file='fes-rew_square_sparse_walls.dat'
    with open(file,'r') as ifile:
        lines=ifile.readlines()
        for i in range(10):
            if ' SET nbins_' in lines[i]:
                size[lines[i].split()[-2].split('_')[1]]=int(lines[i].split()[-1].rstrip())
        print("WALLS")
        print(size)
    
#%%

data=np.loadtxt(file)
print(data[0:-1][0])
cv1_temp=[x[0] for x in data]
cv2_temp=[x[1] for x in data]
free=[x[2] for x in data]

#print(data[3].split()[1])
#for i in range(3):
#    if data[i].split()[1] in data[3].split()[1]:
#        print('ciao')

#for line in data:
#            if float(line.split()[0]) not in cv1:
#                cv1.append(float(line.split()[0]))
                
#            if float(line.split()[1]) not in cv2:
#                cv2.append(float(line.split()[1]))

#cv1=[line.split()[0]  for line in data if line.split()[0] not in cv1]
#cv2=[line.split()[1]  for line in data if line.split()[1] not in cv2]
#cv2=[float(line.split()[0])  for line in lines if line.rstrip() if "#" not in line if line.split()[0] not in cv2 ]
#free=np.array(data[2])
#der_cv1=np.array([float(line.split()[3])  for line in data])
#der_cv2=np.array([float(line.split()[4])  for line in data])
cv1=np.unique(cv1_temp)
cv2=np.unique(cv2_temp)
grid=np.array(free).reshape((len(cv1),len(cv2)))
print(grid.shape)
print(len(cv1),len(cv2),len(free))
print(cv1[1],cv2[0],grid[0][1])
#%%
print(np.where(grid == np.min(grid)))
#print(grid[119][32])
#print(cv1[32],cv2[119])
#%%
fig=plt.figure(figsize=(5,5),dpi=150)
lev=int(round(np.max(grid)/4,0))
plt.contourf(cv1, cv2,grid,lev)
plt.colorbar()
  
# %% PLOT SLICE
minima=np.loadtxt('minima.dat',unpack=True)
slice=np.unique(minima[1])
#slice=[x x-slice[j]]
print(slice)
diff=[[abs(x-k) for x in cv1] for k in slice]
diff2=[[abs(x-k) for x in cv2] for k in slice]
act_slice=[cv1[x.index(min(x))] for x in diff]
act_slice2=[cv2[x.index(min(x))] for x in diff2]
idx=[np.where(cv1 == x) for x in act_slice]
idx2=[np.where(cv2 == x) for x in act_slice2]
print(grid[1][0])
slice_arr=[[grid[x][k] for x in range(len(grid)) ] for k in idx]
slice_arr2=[grid[k][0] for k in idx2]
#equal=range(0.5,2.0,0.3)    
#slice_arr3=[grid[idx2][0]
for i,item in enumerate(slice):
    np.savetxt(f"slice_{item}.dat",list(zip(cv1,slice_arr2[i])),"%.3f")


#print(slice_arr2)
#print(grid[120][32])
#print(min(slice_arr))
#print(min(grid[:][32]))
fig=plt.figure(figsize=(6,5),dpi=150)
#plt.plot(cv2,slice_arr,label="CV1 = {}".format(act_slice))
for i,item in enumerate(slice):
    plt.plot(cv1,slice_arr2[i],label="CV2 = {}".format(act_slice2[i]))
plt.xlim([-2.5,2.5])
plt.ylim([0,150])
plt.xlabel("CV")
plt.ylabel("E / Kj/mol")
plt.legend()
plt.title(os.path.split(os.path.split(file)[-2])[-1])
plt.savefig(f'slice.png',format='png', dpi=150)

# %%
"""
with open(root+'COLVAR_no_explosion') as ifile:
    lines=ifile.readlines()
    time=[float(x.split()[0]) for x in lines[1:]]
    bias=[float(x.split()[-2]) for x in lines[1:]]
print(max(bias),time[bias.index(max(bias))])

"""
# %%
