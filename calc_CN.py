import numpy as np
import glob
import sys

data=np.loadtxt(sys.argv[1],unpack=True)

def get_closest(data,max):
    diffs=[abs(x-float(max)) for x in data[0]]
    limit=diffs.index(min(diffs))
    return(limit+1)

def coord_Number(x,y,dist_min,dist_max):
    val_x=[k for k in x if dist_min<= k <=dist_max]
    val_x_sq=[k**2 for k in val_x]
    val_y=[y[i] for i,k  in enumerate(x) if dist_min <= k <= dist_max]
    prod=[k*val_y[i] for i,k in enumerate(val_x_sq)]
    CN_all=4*np.pi*np.trapz(prod,val_x)
    return(CN_all)


#print(data)
limit=get_closest(data,sys.argv[2])
#print(data[0][:limit])
intg=4*np.pi*np.trapz(data[1][:limit]*data[0][:limit]**2,data[0][:limit])

print(f"Integral up to {sys.argv[2]}: {intg:.3f}")

def calculate_coordination_number(r, g_r, r_min, r_max):
    # Filter the radial distances and RDF values in the specified range
    mask = (r >= r_min) & (r <= r_max)
    r_filtered = r[mask]
    g_r_filtered = g_r[mask]

    # Calculate the coordination number
    dr = np.diff(r_filtered)
    CN = 4 * np.pi * np.sum(g_r_filtered[:-1] * r_filtered[:-1]**2 * dr)
    
    return CN
CN=calculate_coordination_number(data[0],data[1],0,float(sys.argv[2]))
print(f"CN={CN}")