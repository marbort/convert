#%%
import numpy as np
import matplotlib.pyplot as plt
from plot2d_reweight import extract_data
import sys

cv1,cv2,free_grid=extract_data('/home/marco/SHARED/RATIO/WP1/ML/fox/iPrMgCl/SAGA/META/opes/Mg1-O_Mg2-O/40_it25_noexplore/fes-rew_square_sparse.dat')
start=(1.,1.)
end=(3,3)
cvs=np.array([cv1,cv2])


def calc_closest(point:tuple,array):
    diffx=[abs(point[0]-x) for x in array[0]]
    diffy=[abs(point[1]-x) for x in array[1]]
    closest_val=(array[0][np.where(diffx==min(diffx))],array[1][np.where(diffy==min(diffy))])
    closest=(np.where(array[0]==closest_val[0]),np.where(array[1]==closest_val[1]))
    return(closest)

closest_start=calc_closest(start,cvs)
closest_end=calc_closest(end,cvs)


#Initialize auxiliary arrays
distmap=np.ones((len(cv1),len(cv2)),dtype=int)*np.Infinity

map=np.transpose(free_grid)
originmap=np.ones((len(cv1),len(cv2)),dtype=float)*np.nan
visited=np.zeros((len(cv1),len(cv2)),dtype=bool)
finished = False
x,y=closest_start[0][0][0],closest_start[1][0][0]
distmap[x,y]=0
count=0
#print(x,y,map[closest_start[0],closest_start[1]])
#print(np.transpose(free_grid)[closest_start[0],closest_start[1]])
#Loop Dijkstra until reaching the target cell
while not finished:
#while count < 1:
  # move to x+1,y
  if x < closest_end[0][0]-1:
    if distmap[x+1,y]>map[x+1,y]+distmap[x,y] and not visited[x+1,y]:
      distmap[x+1,y]=map[x+1,y]+distmap[x,y]
      originmap[x+1,y]=np.ravel_multi_index([x,y], (closest_end[0][0][0],closest_end[1][0][0]))
      print(originmap[x+1,y])
  # move to x-1,y
  if x>0:
    if distmap[x-1,y]>map[x-1,y]+distmap[x,y] and not visited[x-1,y]:
      distmap[x-1,y]=map[x-1,y]+distmap[x,y]
      originmap[x+1,y]=np.ravel_multi_index([x,y], (closest_end[0][0][0],closest_end[1][0][0]))
      print(originmap[x-1,y])
  # move to x,y+1
  if y < closest_end[1][0]-1:
    if distmap[x,y+1]>map[x,y+1]+distmap[x,y] and not visited[x,y+1]:
      distmap[x,y+1]=map[x,y+1]+distmap[x,y]
      originmap[x+1,y]=np.ravel_multi_index([x,y], (closest_end[0][0][0],closest_end[1][0][0]))
      print(originmap[x,y+1])
  # move to x,y-1
  if y>0:
    if distmap[x,y-1]>map[x,y-1]+distmap[x,y] and not visited[x,y-1]:
      distmap[x,y-1]=map[x,y-1]+distmap[x,y]
      originmap[x+1,y]=np.ravel_multi_index([x,y], (closest_end[0][0][0],closest_end[1][0][0]))
      print(originmap[x,y-1])
  visited[x,y]=True
  dismaptemp=distmap
  dismaptemp[np.where(visited)]=np.Infinity
  #print(dismaptemp[x,y-1])
  # now we find the shortest path so far
  minpost=np.unravel_index(np.argmin(dismaptemp),np.shape(dismaptemp))
  #x=np.min(originmap-)
  x,y=minpost[0],minpost[1]
  if x==closest_end[0][0]-1 and y==closest_end[1][0]-1:
    finished=True
  count=count+1
  #print(x,y)
np.savetxt('/home/marco/SHARED/RATIO/WP1/ML/fox/iPrMgCl/SAGA/META/opes/Mg1-O_Mg2-O/40_it25_noexplore/test.txt',
           originmap,fmt='%.3f')  

#print(np.int(originmap[x,y]))
#Start backtracking to plot the path  
mattemp=map.astype(float)
x,y=closest_end[0][0][0]-1,closest_end[1][0]-1
print(x,y)
path=[]
mattemp[int(x),int(y)]=np.nan

#while x>closest_start[0][0][0] or y>closest_start[1][0][0]:
number=0
while number < 10:
    print("Running {} {}".format(x,y))
    path.append([int(x),int(y)])
    #print(int(originmap[int(x),int(y)]),closest_end[0][0])
    xxyy=np.unravel_index(int(originmap[int(x),int(y)]), (closest_end[0][0][0],closest_end[1][0][0]))
    x,y=xxyy[0],xxyy[1]
    mattemp[int(x),int(y)]=np.nan
    number +=1
path.append([int(x),int(y)])

print(path)
#map[max_val-1,max_val-1]=closest_end

#%%
"""
fig, ax = plt.subplots(figsize=(10,10))

maxnum=np.max(free_grid)
map = free_grid
closest_start=closest(start,[cv1,cv2])
closest_end=closest(end,[cv1,cv2])
min_val, max_val = int(min(cv1)),int(max(cv1))
map[0,0]=closest_start
map[max_val-1,max_val-1]=closest_end

current_cmap = plt.cm.Blues
current_cmap.set_bad(color='red')
ax.matshow(map, cmap=plt.cm.Blues, vmin=0, vmax=maxnum*2)

#Initialize auxiliary arrays
distmap=np.ones((max_val,max_val),dtype=int)*np.Infinity
distmap[0,0]=0
originmap=np.ones((max_val,max_val),dtype=int)*np.nan
visited=np.zeros((max_val,max_val),dtype=bool)
finished = False
x,y=np.int(0),np.int(0)
count=0

#Loop Dijkstra until reaching the target cell
while not finished:
  # move to x+1,y
  if x < max_val-1:
    if distmap[x+1,y]>map[x+1,y]+distmap[x,y] and not visited[x+1,y]:
      distmap[x+1,y]=map[x+1,y]+distmap[x,y]
      originmap[x+1,y]=np.ravel_multi_index([x,y], (max_val,max_val))
  # move to x-1,y
  if x>0:
    if distmap[x-1,y]>map[x-1,y]+distmap[x,y] and not visited[x-1,y]:
      distmap[x-1,y]=map[x-1,y]+distmap[x,y]
      originmap[x-1,y]=np.ravel_multi_index([x,y], (max_val,max_val))
  # move to x,y+1
  if y < max_val-1:
    if distmap[x,y+1]>map[x,y+1]+distmap[x,y] and not visited[x,y+1]:
      distmap[x,y+1]=map[x,y+1]+distmap[x,y]
      originmap[x,y+1]=np.ravel_multi_index([x,y], (max_val,max_val))
  # move to x,y-1
  if y>0:
    if distmap[x,y-1]>map[x,y-1]+distmap[x,y] and not visited[x,y-1]:
      distmap[x,y-1]=map[x,y-1]+distmap[x,y]
      originmap[x,y-1]=np.ravel_multi_index([x,y], (max_val,max_val))

  visited[x,y]=True
  dismaptemp=distmap
  dismaptemp[np.where(visited)]=np.Infinity
  # now we find the shortest path so far
  minpost=np.unravel_index(np.argmin(dismaptemp),np.shape(dismaptemp))
  x,y=minpost[0],minpost[1]
  if x==max_val-1 and y==max_val-1:
    finished=True
  count=count+1

#Start backtracking to plot the path  
mattemp=map.astype(float)
x,y=max_val-1,max_val-1
path=[]
mattemp[np.int(x),np.int(y)]=np.nan

while x>0.0 or y>0.0:
  path.append([np.int(x),np.int(y)])
  xxyy=np.unravel_index(np.int(originmap[np.int(x),np.int(y)]), (max_val,max_val))
  x,y=xxyy[0],xxyy[1]
  mattemp[np.int(x),np.int(y)]=np.nan
path.append([np.int(x),np.int(y)])
print(path)

#Output and visualization of the path
current_cmap = plt.cm.Blues
current_cmap.set_bad(color='red')
fig, ax = plt.subplots(figsize=(8,8))
ax.matshow(mattemp,cmap=plt.cm.Blues, vmin=0, vmax=20)
for i in range(max_val):
    for j in range(max_val):
      c = map[j,i]
      ax.text(i, j, str(c), va='center', ha='center')
plt.savefig('minpath.png',format='png')

print('The path length is: '+np.str(distmap[max_val-1,max_val-1]))
print('The dump/mean path should have been: '+np.str(maxnum*max_val))

"""

