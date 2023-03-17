#%%
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

root='/home/marco/SHARED/RATIO/WP1/ML/fox/iPrMgCl/SAGA/Meta_test/NEW'
slice='1.509412933'
cv=2
size={}
with open(os.path.join(root,'fes.dat'),'r') as ifile:
    lines=ifile.readlines()
    for i in range(10):
        if ' SET nbins_' in lines[i]:
            size[lines[i].split()[-2].split('_')[1]]=int(lines[i].split()[-1].rstrip())
    print(size)
#%%
cv1=[]
cv2=[]
data=[x for x in lines if "#" not in x if x.rstrip()]
print(data[3].split()[1])
for i in range(3):
    if data[i].split()[1] in data[3].split()[1]:
        print('ciao')

for line in data:
            if float(line.split()[0]) not in cv1:
                cv1.append(float(line.split()[0]))
                
            if float(line.split()[1]) not in cv2:
                cv2.append(float(line.split()[1]))

#cv1=[line.split()[0]  for line in data if line.split()[0] not in cv1]
#cv2=[line.split()[1]  for line in data if line.split()[1] not in cv2]
#cv2=[float(line.split()[0])  for line in lines if line.rstrip() if "#" not in line if line.split()[0] not in cv2 ]
free=np.array([float(line.split()[2])  for line in data])
der_cv1=np.array([float(line.split()[3])  for line in data])
der_cv2=np.array([float(line.split()[4])  for line in data])

grid=free.reshape((len(cv1),len(cv2)))
print(grid.shape)
print(len(cv1),len(cv2),len(free))
print(cv1[1],cv2[0],grid[0][1])
#%%
print(np.where(grid == np.min(grid)))
print(grid[119][32])
print(cv1[32],cv2[119])
#%%
fig=plt.figure(figsize=(5,5),dpi=150)
lev=int(round(np.max(grid)/4,0))
plt.contourf(cv1, cv2,grid,lev)
plt.colorbar()
  
# %% PLOT SLICE
slice=1.00000000
diff=[abs(x-slice) for x in cv1]
act_slice=cv1[diff.index(min(diff))]
idx=cv1.index(act_slice)
slice_arr=[grid[x][idx] for x in range(len(grid)) ]
print(idx)
print(grid[120][32])
print(min(slice_arr))
#print(min(grid[:][32]))
plt.plot(cv2,slice_arr)
plt.xlabel("CV2")
plt.ylabel("E / Kj/mol")
plt.title("CV1 = {}".format(act_slice))

# %%
with open(root+'COLVAR_no_explosion') as ifile:
    lines=ifile.readlines()
    time=[float(x.split()[0]) for x in lines[1:]]
    bias=[float(x.split()[-2]) for x in lines[1:]]
print(max(bias),time[bias.index(max(bias))])


# %%
