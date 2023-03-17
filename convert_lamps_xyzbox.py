#%%
import dpdata as dp
import glob
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import regex as re
from natsort import natsorted, ns
import json


root='/home/marco/SHARED/RATIO/WP1/ML/fox/iPrMgCl/SAGA/'
sys="iPrMgCl-dimer"
natoms=2572

#string=root+'Iteration*/'+sys+'_*/lmp*'
paths=glob.glob(root+'MD_test/lmp*')
paths_sort=natsorted(paths, key=lambda y: y.lower())
#dir_path = os.path.dirname(os.path.realpath(c[0]))
print(paths_sort)
#print(re.findall('Iteration_\d*',dir_path))

def create_plumed_input(root,name):
    with open(os.path.join(root,name+'_input.json'),'r') as ifile:
        input_dict=json.load(ifile)
    with open(os.path.join(root,name+'.dat'),'w') as pf:
        tset=(os.path.dirname(os.path.realpath(paths_sort[0])) + "/" + "conf.data")
        #print(tset,args.meta,previous_iteration)
        data=dp.System(tset,'lammps/lmp')
        groups=[np.where(data['atom_types']==data['atom_names'].index(input_dict['plumed']['groups'][x][0])) for x in input_dict['plumed']['groups']]
        for i,item in enumerate(input_dict['plumed']['groups']):
            if input_dict['plumed']['groups'][item][1]=="all":
                sel_group=groups[i][0]
                sel_group=[x+1 for x in sel_group]
                sel_group_str=np.char.mod('%0d', sel_group)

                pf.write(item+": GROUP ATOMS="+','.join(sel_group_str)+'\n')
            else:
                sel_group=[groups[i][0][input_dict['plumed']['groups'][item][1][x]] for x in range(len(input_dict['plumed']['groups'][item][1]))]
                sel_group=[x+1 for x in sel_group]
                sel_group_str=np.char.mod('%0d', sel_group)
                pf.write(item+": GROUP ATOMS="+','.join(sel_group_str)+'\n')
            pf.write("\n\n")
        for i in input_dict['plumed']['cvs']:
                pf.write(i+'\n')
        pf.write("\n\n")
        
        pf.write(input_dict['plumed']['print']+'\n\n')


#%%
RT="/home/marco/SHARED/RATIO/WP1/ML/fox/iPrMgCl/SAGA/MD_test"
plmd="analysis_iPrMgCl-dimer_Mg-O"
create_plumed_input(RT,plmd)
#%%
for path in paths_sort:
    dump_freq=0
    dir_path = os.path.dirname(os.path.realpath(path))
    if os.path.exists(dir_path+'/colvar_number2'):
        print("\nSkipping {}.\n".format(dir_path))
    else:
    
        data=dp.System(path,'lammps/dump')
        box=[[str(x[i][i]) for i in range(3)] for x in data['cells']]
        data.to('xyz',dir_path+'/trj.xyz')

        with open(dir_path+'/trj.xyz','r') as ifile:
            lines=ifile.readlines()
            k=0
            for i,line in enumerate(lines):
                if line =='\n':
                    lines[i]=" ".join(box[k])+"\n"
                    k=k+1

        with open(dir_path+'/trj_box.xyz','w') as ofile:
            for line in lines:
                ofile.write(line)
    
        with open(dir_path+'/start.lmp') as ifile:
            lines=ifile.readlines()
            for line in lines:
                if 'dump_freq equal' in line:
                    dump_freq=line.split()[-1]
        
        
        subprocess.call(['plumed','driver','--plumed',root+plmd+'.dat','--ixyz',dir_path+'/trj_box.xyz',\
                        '--timestep', '0.0001', '--trajectory-stride', dump_freq, '--length-units', 'A'])
        subprocess.call(['cp','colvar_both',dir_path])
        print("Finished "+path)
print("\n\nEND!!\n")
#%%

CVs={}
fields=[]
for path in paths_sort:
    dir_path = os.path.dirname(os.path.realpath(path))
    it=re.findall('Iteration_\d*',dir_path)[0]
    CVs[it]={}
print(CVs)
for path in paths_sort:
    dir_path = os.path.dirname(os.path.realpath(path))
    it=re.findall('Iteration_\d*',dir_path)[0]
    env=dir_path.split('/')[-1].split('_')[-2]+" "+dir_path.split('/')[-1].split('_')[-1]
    print(env)
    with open(dir_path+'/colvar') as ifile:
        lines=ifile.readlines()
        fields=[lines[0].split()[i] for i in range(2,len(lines[0].split()))]
        cond={env:{x:[float(line.split()[i]) for line in lines[1:]] for i,x in enumerate(fields)}}
        CVs[it].update(cond)

#%%
for i in CVs:
    print(i,len(CVs[i]))
    
#%%      CV1
fig=plt.figure(figsize=(20,15),dpi=150)
cols=4
rows=(len(CVs)//4+1)
print(CVs['Iteration_1'])
for i,item in enumerate(CVs):
    lgnd=[]
    plt.subplot(rows,cols,i+1)
    for env in CVs[item]:
        cond="{}".format(env)
        plt.plot(CVs[item][env][fields[0]],CVs[item][env][fields[1]])
        lgnd.append(env)
    plt.legend(lgnd)
    plt.title(i)


#%%      CV2
fig=plt.figure(figsize=(20,15),dpi=150)
cols=4
rows=(len(CVs)//4+1)
print(CVs['Iteration_1'])
for i,item in enumerate(CVs):
    lgnd=[]
    plt.subplot(rows,cols,i+1)
    item="Iteration_{}".format(i+1)
    for env in CVs[item]:
        cond="{}".format(env)
        plt.plot(CVs[item][env][fields[0]],CVs[item][env][fields[2]])
        lgnd.append(env)
    plt.legend(lgnd)
    plt.title(item)

    
    #plumed driver --plumed analysis.dat --ixyz trj_box --timestep 0.0001 --trajectory-stride 1000 --length-units A
            
    

# %%BOTH per env
#fig=plt.figure(figsize=(20,15),dpi=150)
cols=1
rows=(len(CVs)//cols+1)
#print(CVs['Iteration_1'])
interesting_env=["1.0bar 200K","1.0bar 300K"]
for j in interesting_env:
    plt.figure(figsize=(17,22),dpi=150)
    plt.suptitle(j,y=1.00)
    limx=[0,25]
    limy=[1,3]
    for i,item in enumerate(CVs):
        lgnd=[]
        plt.subplot(rows,cols,i+1)
        for k in range(1,len(fields)):
                plt.plot(CVs[item][j][fields[0]],CVs[item][j][fields[k]])
                lgnd.append(fields[k])
        plt.legend(lgnd)
        plt.title(item)
        plt.xlim(limx)
        plt.ylim(limy)
    plt.tight_layout()
# %%
