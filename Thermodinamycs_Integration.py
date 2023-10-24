import os
import numpy as np
import glob
import matplotlib.pyplot as plt



def extract_Langrange_Multipliers(files):
    data={}
    for file in files: 
            with open(file,'r') as ifile:
                    dir=os.path.dirname(file)
                    lines=ifile.readlines()
                    forces=[float(x.split()[-1]) for x in lines[1:] if "Shake" in x]
                    avg_force=np.mean(forces)
                    std_err=np.std(forces)/np.sqrt(len(forces))
                    cv=float(dir.split('_')[-1])
                    data[dir]={'cv':cv,'frcs':forces,'avg_frc':avg_force,'std_err':std_err}
                    
    return(data)

def plot_Langrange_Multipliers(data,wd,hg):
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    maxim=0
    minim=0
    for cv in data:
        if abs(max(data[cv]['frcs'])) > abs(maxim):
            maxim=max(data[cv]['frcs'])
        if abs(min(data[cv]['frcs'])) > abs(minim):
            minim=min(data[cv]['frcs'])
            
    for i,cv in enumerate(data):
        print(cv)
        plt.subplot(len(data)//3+1,3,i+1)
        plt.plot(list(range(1,len(data[cv]['frcs'])+1)),data[cv]['frcs'])
        plt.title(cv)
        plt.ylim([minim-abs(minim*0.1),maxim+abs(maxim*0.1)])
    plt.tight_layout()
    plt.savefig('LagrMult.png')
    
   
    return()



def integrate_force(data):
    FES=[]
    cv=[(data[x]['cv'],data[x]['avg_frc'],data[x]['std_err']) for x in data]
    cv.sort()
    print(cv)
    for i in range(len(cv)):
        FES.append(-np.trapz([cv[0][1],cv[i][1]],[cv[0][0],cv[i][0]]))
    with open('FES.dat','w') as ofile:
        ofile.write("{:10s} {:10s} {:10s} {:10s} {:10s}".format("#VAR","AVG_FRC","FRC_ERR","FES / Ha","FES / Kjmol-1\n"))
        for i,var in enumerate(cv):
            ofile.write("{:10.5f} {:10.5f} {:10.5f} {:10.5f} {:10.5f}\n".format(var[0]*2,var[1],var[2],FES[i],FES[i]*2625.5))
    return(cv,FES)     
def plot_avg_force_TI(data,FES,cv):
    plt.subplot(2,1,1)
    plt.errorbar([data[i]['cv'] for i in data],[data[i]['avg_frc'] for i in data],yerr=[data[i]['std_err'] for i in data])
    plt.subplot(2,1,2)
    plt.plot([x[0] for x in cv],FES,'-o')
    plt.savefig('TI.png')
    
    

files=glob.glob("*/*.LagrangeMultLog")
data=extract_Langrange_Multipliers(files)
plot_Langrange_Multipliers(data,15,10)
cv,FES=integrate_force(data)
#print(data)
print("FES {}".format(FES))
plot_avg_force_TI(data,[x* 2625.5 for x in FES],cv)
