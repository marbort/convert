import matplotlib.pyplot as plt
import sys

forcex=[]
forcey=[]
forcez=[]

forcex_set=[]
forcey_set=[]
forcez_set=[]


with open(sys.argv[1],'r') as ifile:
    lines=ifile.readlines()
    for line in lines:
        if line[0].isdigit():
            continue
        if line[0:7] == 'Lattice':
            continue
        else:
            forcex.append(float(line.split()[8]))
            forcey.append(float(line.split()[9]))
            forcez.append(float(line.split()[10]))
            
            forcex_set.append(float(line.split()[4]))
            forcey_set.append(float(line.split()[5]))
            forcez_set.append(float(line.split()[6]))

force_set = [( forcex_set[i]**2 + forcey_set[i]**2 + forcez_set[i]**2 )**0.5 for i in range(len(forcex_set)) ]
force = [( forcex[i]**2 + forcey[i]**2 + forcez[i]**2 )**0.5 for i in range(len(forcex)) ]

fig, ax = plt.subplots(2, 2, figsize=(10, 10), dpi=150)
size=1  
ax[0, 0].scatter(forcex_set,forcex, label='Force X', color='blue',s=size)
ax[0, 0].set_title('Ref vs Predicted Z')
ax[0, 0].set_xlabel('Ref')
ax[0, 0].set_ylabel('Predicted')
ax[0, 0].legend()
ax[0, 1].scatter(forcey_set,forcey, label='Set Y Force', color='orange',s=size)
ax[0, 1].set_title('Ref vs Predicted Z')
ax[0, 1].set_xlabel('Ref')
ax[0, 1].set_ylabel('Predicted')
ax[0, 1].legend()
ax[1, 0].scatter(forcez_set,forcez, label='Set Z Force', color ='green',s=size)
ax[1, 0].set_title('Ref vs Predicted Z')
ax[1, 0].set_xlabel('Ref')
ax[1, 0].set_ylabel('Predicted')
ax[1, 0].legend()
ax[1, 1].scatter(force_set,force, label='Total Force', color='red',s=size)
ax[1, 1].set_title('Ref vs Predicted Total Force')
ax[1, 1].set_xlabel('Ref')
ax[1, 1].set_ylabel('Predicted')
ax[1, 1].legend()

plt.tight_layout()
plt.savefig('force_comparison.png')
plt.close(fig)
    
            
