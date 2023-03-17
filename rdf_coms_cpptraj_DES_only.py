#%%
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import os
from MDAnalysis.analysis import *
from MDAnalysis.analysis.bat import BAT
from MDAnalysis.analysis.density import DensityAnalysis
from MDAnalysis.analysis.lineardensity import LinearDensity
import glob
from scipy.signal import argrelextrema
from scipy.signal import find_peaks

#%%
def coord_Number(x,y,dist_min,dist_max,mols):
    val_x=[k for k in x if dist_min<= k <=dist_max]
    val_x_sq=[k**2 for k in val_x]
    val_y=[y[i] for i,k  in enumerate(x) if dist_min <= k <= dist_max]
    prod=[k*val_y[i] for i,k in enumerate(val_x_sq)]
    CN_all=4*np.pi*np.trapz(prod,val_x)
    CN_one=CN_all/mols
    return(CN_one)
def coord_number2(x,y,dist_min,dist_max):
    val_x=[k for k in x if dist_min<= k <=dist_max]
    
    val_y=[y[i] for i,k  in enumerate(x) if dist_min <= k <= dist_max]
    
    CN_all=np.trapz(val_y,val_x)
    return(CN_all)
#%%
with open(root+'cpptraj.atomic.out') as ifile:
    lines=ifile.readlines()
    names=[]
    Dd=[]
    numb=[]

    for i,line in enumerate(lines):
        if "Calculating RDF" in line:
            names.append(line.split()[7]+line.split()[-1])
        if "in common" in line:
           numb.append(max(int(line.split()[6].rstrip(',')),int(line.split()[11].rstrip(','))))
        if "Average density" in line:
            Dd.append(float(line.split()[3]))
        
D={}
#print(names,Dd)
for i,x in enumerate(ords):
    D[x]={'name':names[i],'density':round(Dd[i],6)}

print(D['0024']['density'])
#%%
print(RDFs['0001_CHL_N1_CHL_N1']['CN'],RDFs['0001_CHL_N1_CHL_N1']['r_min'])
print(RDFs['0024_GCL_C2_GCL_C2']['CN'],RDFs['0024_GCL_C2_GCL_C2']['r_min'])
print(RDFs['0025_GCL_C2_GCL_C1C3']['CN'],RDFs['0025_GCL_C2_GCL_C1C3']['r_min'])
print(RDFs['0036_Clm_Clm']['CN'],RDFs['0036_Clm_Clm']['r_min'])
# %%
resnumb={"CHL":400,"Clm":400,"GCL":800}
root="/home/marco/SHARED/RATIO/WP4/MD/MOD-FRC/BIG/DES/"
id="_rdf_vol.dat"
coms=glob.glob(root+"*"+id)
coms.sort()
RDFs={}
names=[]
RESs=[]
mols_n=[]
ords=[]
r_min_lit=[6.4,8.1,7.2,5.6,6.8,7.2,4.5,6.5,5.2,3.4,4.2,4.3,3.5,7.8,8.3,3.4,4.2,4.1,4.3,2.5,4.7,3.4,7.5,7.5,4.3,7.2,7.4,7.3,3.2,3.1,2.4,3.05,2.3,3.4,10.3]
for i in coms:
    with open(i,'r') as ifile:
        name=os.path.basename(i).replace(id,"")
        names.append(name)
        ords.append(name.split('_')[0])
        lines=ifile.readlines()
        res=name.split('_')[1]
        nmols=len(name.split('_')[2].split(','))
        x=[]
        y=[]
        y_min=[]
        size=len(lines[0].split()[1].split('=')[0].split(','))
        for line in lines[1:]:
                x.append(float(line.split()[0]))
                y.append(float(line.split()[1]))
        y_max_index=find_peaks(y,width=1)
        x_max=x[y_max_index[0][0]]
        y_min=find_peaks([-x for x in y],width=1,prominence=0.05)
        x_min=x[y_min[0][0]]
       
            
        
    RDFs[name.split('_')[0]]={'name':'_'.join(name.split('_')[1:]),'size':size,'RDF':[x,y],'r_max':x_max,'r_min':x_min,'minima':y_min,'maxima':y_max_index,'#_part':nmols*resnumb[res]}


for i,x in enumerate(RDFs):
    RDFs[x]['CN']=coord_Number(RDFs[x]['RDF'][0],RDFs[x]['RDF'][1],0,RDFs[x]['r_min'],RDFs[x]['#_part'])*D[x]['density']/RDFs[x]['size']
    RDFs[x]['CN_rmin_lit']=coord_Number(RDFs[x]['RDF'][0],RDFs[x]['RDF'][1],0,r_min_lit[i],RDFs[x]['#_part'])*D[x]['density']/RDFs[x]['size']
#print(RDFs['GCL_C2_GCL_C2'])
with open (root+'results.dat','w') as ofile:
    ofile.write("{:35s} {:>15s} {:>15s} {:>15s} {:>15s}\n".format('RDF','r_max','r_min','CN','CN_r_min_lit'))
    for i in RDFs:
        ofile.write("{:35s} {:15.3f} {:15.3f} {:15.3f} {:15.3f}\n".format(RDFs[i]['name'],RDFs[i]['r_max'],RDFs[i]['r_min'],RDFs[i]['CN'],RDFs[i]['CN_rmin_lit']))
#%%
print(RDFs['0001']['name'],RDFs['0001']['CN'],RDFs['0001']['size'],RDFs['0001']['r_min'])


#%%
print(RDFs['0001_CHL_N1_CHL_N1']['CN'],RDFs['0001_CHL_N1_CHL_N1']['r_min'])
print(RDFs['0024_GCL_C2_GCL_C2']['CN'],RDFs['0024_GCL_C2_GCL_C2']['r_min'])
print(RDFs['0025_GCL_C2_GCL_C1C3']['CN'],RDFs['0025_GCL_C2_GCL_C1C3']['r_min'])
print(RDFs['0036_Clm_Clm']['CN'],RDFs['0036_Clm_Clm']['r_min'])



coord_number(RDFs['0001_CHL_N1_CHL_N1']['RDF'][0],RDFs['0001_CHL_N1_CHL_N1']['RDF'][1],0,6.4,400)
#%%
pipo=coord_Number(RDFs['0024_GCL_C2_GCL_C2']['RDF'][0],RDFs['0024_GCL_C2_GCL_C2']['RDF'][1],0,7.5,800)
print(D['0024']['density']*pipo)
plt.plot(RDFs['0024_GCL_C2_GCL_C2']['RDF'][0],RDFs['0024_GCL_C2_GCL_C2']['RDF'][1])
plt.scatter(RDFs['0024_GCL_C2_GCL_C2']['r_max'],RDFs['0024_GCL_C2_GCL_C2']['RDF'][1][RDFs['0024_GCL_C2_GCL_C2']['maxima'][0][0]])
plt.scatter(RDFs['0024_GCL_C2_GCL_C2']['r_min'],RDFs['0024_GCL_C2_GCL_C2']['RDF'][1][RDFs['0024_GCL_C2_GCL_C2']['minima'][0][0]])
# %% PLOT RDFS
# %% PLOT RDFS
#fig=plt.figure(dpi=150)

#print(np.array(RDFs['CHL_CHL']['RDF'][1]).shape)
cell='0036_Clm_Clm'
#minima=find_peaks([-x for x in RDFs[cell]['RDF'][1]],prominence=0.5)
#print(minima[0])
#print(RDFs[cell]['RDF'][0][minima[0][0]])
plt.plot(RDFs[cell]['RDF'][0],RDFs[cell]['RDF'][1])
#peaks=[[RDFs[cell]['RDF'][0][x],RDFs[cell]['RDF'][1][x]] for x in minima[0]]
plt.scatter(peaks[0][0],peaks[0][1])
#plt.xlim([1.5,2.5])
plt.show()
print(RDFs[cell]['r_max'])
print(RDFs[cell]['r_min'])
#%%
cols=5
rows=len(RDFs)//cols+1

limx=[1,11]
fig=plt.figure(figsize=(30,25),dpi=150)
for i,item in enumerate(RDFs):
    plt.subplot(rows,cols,i+1)
    plt.plot(RDFs[item]['RDF'][0],RDFs[item]['RDF'][1])
    plt.scatter(RDFs[item]['r_max'],RDFs[item]['RDF'][1][RDFs[item]['maxima'][0][0]],s=10)
    plt.scatter(RDFs[item]['r_min'],RDFs[item]['RDF'][1][RDFs[item]['minima'][0][0]],s=10)
    plt.title(RDFs[item]['name'])
    plt.xlim(limx)
plt.tight_layout()
plt.show()

    
    
    

# %%

with open(root+'cpptraj.out') as ifile:
    lines=ifile.readlines()
    names=[]
    Dd=[]
    numb=[]
    for i,line in enumerate(lines):
        if "Calculating RDF" in line:
            names.append(line.split()[7]+line.split()[-1])
        if "in common" in line:
           numb.append(max(int(line.split()[6].rstrip(',')),int(line.split()[11].rstrip(','))))
        if "Average density" in line:
            Dd.append(float(line.split()[3]))
        
D2={}

id2="_rdf.dat"
coms2=glob.glob(root+"*com*"+id2)
RDFs2={}
names2=[]
ords2=[]
print(coms2)
for i in coms2:
    with open(i,'r') as ifile:
        name=os.path.basename(i).replace(id2,"")
        names2.append(name)
        lines=ifile.readlines()
        res=name.split('_')[2]
        ords2.append(name.split('_')[0])
        x=[]
        y=[]
        y_min=[]
        for line in lines[1:]:
                x.append(float(line.split()[0]))
                y.append(float(line.split()[1]))
        y_max_index=find_peaks(y,width=1)
        x_max=x[y_max_index[0][0]]
        y_min=find_peaks([-x for x in y],width=2)
        x_min=x[y_min[0][0]]
    

    
           
        
    RDFs2[name.split('_')[0]]={'name':'_'.join(name.split('_')[1:]),'RDF':[x,y],'r_max':x_max,'r_min':x_min,'minima':y_min,'maxima':y_max_index,'#_part':resnumb[res]}
for i,x in enumerate(ords2):
    D2[x]={'name':names[i],'density':round(Dd[i],6)}

for i in RDFs2:
    RDFs2[i]['CN']=coord_Number(RDFs2[i]['RDF'][0],RDFs2[i]['RDF'][1],0,RDFs2[i]['r_min'],RDFs2[i]['#_part'])*D2[i]['density']
#print(RDFs['GCL_C2_GCL_C2'])

with open (root+'results_COM.dat','w') as ofile:
    ofile.write("{:35s} {:>15s} {:>15s} {:>15s}\n".format('RDF','r_max','r_min','CN'))
    for i in RDFs2:
        ofile.write("{:35s} {:15.3f} {:15.3f} {:15.3f}\n".format(RDFs2[i]['name'],RDFs2[i]['r_max'],RDFs2[i]['r_min'],RDFs2[i]['CN']))

# %%
cell2='com_CHL_Clm'
plt.plot(RDFs2[cell2]['RDF'][0],RDFs2[cell2]['RDF'][1])
#plt.xlim([6,7])
plt.show()
print(RDFs2[cell2]['r_max'])
print(RDFs2[cell2]['r_min'])

cols=3
rows=len(RDFs2)//cols+1
limy2=[0,8]
limx2=[0,14]
fig=plt.figure(figsize=(20,15),dpi=150)
for i,item in enumerate(RDFs2):
    plt.subplot(rows,cols,i+1)
    plt.plot(RDFs2[item]['RDF'][0],RDFs2[item]['RDF'][1])
    plt.title(item)
    plt.xlim(limx2)
    plt.ylim(limy2)
xb_int=[x for x in RDFs2['com_CHL_Clm']['RDF'][0] if x < 6.2]
yb_int=[RDFs2['com_CHL_Clm']['RDF'][1][i] for i,x in enumerate(RDFs2['com_CHL_Clm']['RDF'][0]) if x < 6.2]
print(yb_int)
b=np.trapz(yb_int,xb_int)
print("{} pippo".format(b))
for i in RDFs2:
    print(i, RDFs2[i]['CN'])

# %%
fig=plt.figure(figsize=(10,10),dpi=150)
clrs=['blue','orange','green','black','red','purple']
legd=["Choline-Choline","Choline-chloride","Choline-glycerol",'Chloride-Chloride',"Glycerol-chloride","Glycerol-glycerol",'Chloride-Chloride']
for i,item in enumerate(RDFs2):
    plt.plot(RDFs2[item]['RDF'][0], RDFs2[item]['RDF'][1],color=clrs[i])
plt.legend(legd)    
for i,item in enumerate(RDFs2):
    plt.scatter(RDFs2[item]['r_max'],RDFs2[item]['RDF'][1][RDFs2[item]['maxima'][0][0]],color=clrs[i])
    plt.scatter(RDFs2[item]['r_min'],RDFs2[item]['RDF'][1][RDFs2[item]['minima'][0][0]],color=clrs[i])

plt.title("COM-RDFs")
for i in RDFs2:
    print(RDFs2[i]['name'], RDFs2[i]['CN'])
#plt.xlim(limx2)
plt.ylim([0,3.5])
# %%
id3="_rdf.dat"
coms3=glob.glob(root+"dens*"+id3)
RDFs3={}
names3=[]
print(coms3)
for i in coms3:
    with open(i,'r') as ifile:
        name=os.path.basename(i).replace(id3,"")
        names3.append(name)
        lines=ifile.readlines()
        x=[]
        y=[]
        y_min=[]
        for line in lines:
            if "#" in line:
                continue
            else:
                x.append(float(line.split()[0]))
                y.append(float(line.split()[1]))
        y_max_index=find_peaks(y,width=1)
        x_max=x[y_max_index[0][0]]
        y_min=find_peaks([-x for x in y],width=2)
        x_min=x[y_min[0][0]]
            
        
    RDFs3[name]={'RDF':[x,y],'r_max':x_max,'r_min':x_min,'minima':y_min,'maxima':y_max_index}


for i in RDFs3:
    RDFs3[i]['CN']=coord_number2(RDFs3[i]['RDF'][0],RDFs3[i]['RDF'][1],0,RDFs3[i]['r_min'])
#print(RDFs['GCL_C2_GCL_C2'])

with open (root+'results_COM.dat','w') as ofile:
    ofile.write("{:35s} {:>15s} {:>15s} {:>15s}\n".format('RDF','r_max','r_min','CN'))
    for i in RDFs3:
        ofile.write("{:35s} {:15.3f} {:15.3f} {:15.3f}\n".format(i,RDFs3[i]['r_max'],RDFs3[i]['r_min'],RDFs3[i]['CN']))

# %%
cell2='dens_CHL_Clm'
plt.plot(RDFs3[cell2]['RDF'][0],RDFs3[cell2]['RDF'][1])
#plt.xlim([6,7])
plt.show()
print(RDFs3[cell2]['r_max'])
print(RDFs3[cell2]['r_min'])

cols=3
rows=len(RDFs3)//cols+1
limy3=[0,3.5]
limx3=[0,14]
fig=plt.figure(figsize=(15,10),dpi=150)
for i,item in enumerate(RDFs3):
    plt.subplot(rows,cols,i+1)
    plt.plot(RDFs3[item]['RDF'][0],RDFs3[item]['RDF'][1])
    
    plt.title(item)
    plt.xlim(limx2)
    plt.ylim(limy2)
# %%
fig=plt.figure(figsize=(5,5),dpi=150)
clrs=['blue','orange','green','red','purple']
#legd=["Choline-Choline","Choline-chloride","Choline-glycerol","Glycerol-chloride","Glycerol-glycerol"]
legd=[]
for i,item in enumerate(RDFs3):
    plt.plot(RDFs3[item]['RDF'][0],RDFs3[item]['RDF'][1],color=clrs[i])
    legd.append(item)
plt.legend(legd)
for i,item in enumerate(RDFs3):  
    plt.scatter(RDFs3[item]['r_max'],RDFs3[item]['RDF'][1][RDFs3[item]['maxima'][0][0]],color=clrs[i])
    plt.scatter(RDFs3[item]['r_min'],RDFs3[item]['RDF'][1][RDFs3[item]['minima'][0][0]],color=clrs[i])
plt.title("COM-RDFs")
plt.xlim(limx3)
plt.ylim(limy3)
x_int=[x for x in RDFs3['dens_CHL_Clm']['RDF'][0] if x < 6.2]
y_int=[RDFs3['dens_CHL_Clm']['RDF'][1][i] for i,x in enumerate(RDFs3['dens_CHL_Clm']['RDF'][0]) if x < 6.2]
a=np.trapz(y_int,x_int)
print("{} ciao".format(a))

for i in RDFs3:
    print(i, RDFs3[i]['CN'])

#%%
with open(root+'dens_DES.dat','r') as ifile:
    lines=ifile.readlines()
    z=[float(x.split()[0]) for x in lines if "#" not in x]
    dens=[float(x.split()[1]) for x in lines if "#" not in x]
dens_des=[dens[i] for i,x in enumerate(z) if 10<x<90 ]
avg_dens=np.mean(dens_des)
xy=4054.28
amu_to_g=1.66054e-24
A3_to_mL=1e-24
avg_dens_human=(avg_dens/xy) * amu_to_g /A3_to_mL
print(avg_dens_human)
#plt.plot(z,dens) 
       

# %%
# %%
id4="_rdf.dat"
coms4=glob.glob(root+"raw*"+id3)
RDFs4={}
names4=[]
print(coms4)
for i in coms4:
    with open(i,'r') as ifile:
        name=os.path.basename(i).replace(id4,"")
        names4.append(name)
        lines=ifile.readlines()
        x=[]
        y=[]
        y_min=[]
        for line in lines:
            if "#" in line:
                continue
            else:
                x.append(float(line.split()[0]))
                y.append(float(line.split()[1]))
        y_max_index=find_peaks(y,width=1)
        x_max=x[y_max_index[0][0]]
        y_min=find_peaks([-x for x in y],width=2)
        x_min=x[y_min[0][0]]
            
        
    RDFs4[name]={'RDF':[x,y],'r_max':x_max,'r_min':x_min}


for i in RDFs4:
    RDFs4[i]['CN']=coord_number2(RDFs4[i]['RDF'][0],RDFs4[i]['RDF'][1],0,RDFs4[i]['r_min'])
#print(RDFs['GCL_C2_GCL_C2'])

with open (root+'results_COM.dat','w') as ofile:
    ofile.write("{:35s} {:>15s} {:>15s} {:>15s}\n".format('RDF','r_max','r_min','CN'))
    for i in RDFs4:
        ofile.write("{:35s} {:15.3f} {:15.3f} {:15.3f}\n".format(i,RDFs4[i]['r_max'],RDFs4[i]['r_min'],RDFs4[i]['CN']))

# %%
cell2='dens_CHL_Clm'
plt.plot(RDFs4[cell2]['RDF'][0],RDFs4[cell2]['RDF'][1])
#plt.xlim([6,7])
plt.show()
print(RDFs4[cell2]['r_max'])
print(RDFs4[cell2]['r_min'])

cols=3
rows=len(RDFs4)//cols+1
limy3=[0,3.5]
limx3=[0,14]
fig=plt.figure(figsize=(15,10),dpi=150)
for i,item in enumerate(RDFs4):
    plt.subplot(rows,cols,i+1)
    plt.plot(RDFs4[item]['RDF'][0],RDFs4[item]['RDF'][1])
    plt.title(item)
    plt.xlim(limx2)
    plt.ylim(limy2)
# %%
fig=plt.figure(figsize=(5,5),dpi=150)
clrs=['blue','orange','green','red','purple']
legd=["Choline-Choline","Choline-chloride","Choline-glycerol","Glycerol-chloride","Glycerol-glycerol"]
for i,item in enumerate(RDFs4):
    plt.plot(RDFs4[item]['RDF'][0],RDFs4[item]['RDF'][1],color=clrs[i])

plt.legend(legd)
plt.title("COM-RDFs")
#plt.xlim(limx3)
#plt.ylim(limy3)
x_int=[x for x in RDFs4['dens_CHL_Clm']['RDF'][0] if x < 6.2]
y_int=[RDFs4['dens_CHL_Clm']['RDF'][1][i] for i,x in enumerate(RDFs4['dens_CHL_Clm']['RDF'][0]) if x < 6.2]
a=np.trapz(y_int,x_int)
print("{} ciao".format(a))

for i in RDFs4:
    print(i, RDFs4[i]['CN'])

#%%

#plt.plot(z,dens) 
with open(root+'rdf.xvg') as ifile:
    lines=ifile.readlines()
    x=[float(x.split()[0]) for x in lines if "#" not in x if "@" not in x] 
    y=[float(x.split()[1]) for x in lines if "#" not in x if "@" not in x]
plt.plot(x,y)
coord_number2(x,y,0,8)
       

# %%

# %%
