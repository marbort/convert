#%%
import matplotlib.pyplot as plt
import sys, os
import glob
import matplotlib.colors as mcolors

root="/Volumes/My Passport/RATIO/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/"
files=glob.glob(root+"*/umbrella_30_200/profile.xvg")
histo=glob.glob(root+"*/umbrella_30_200/histo.xvg")
dens=glob.glob(root+"ACT/dens*.xvg")
data=[]
data_hist=[]
chars=["#","@"]
names=[]
rho=[]
labs={"CHL":"CholCl","GCL":"Glycerol","THF":"THF","IPM":"iPr2Mg",
"IMC":"iPrMgCl","MG":"MgCl-","ACT":"Acetone"}
for file in files:
    with open(file,'r') as ifile:
        lines=ifile.readlines()
        #x=[]
        #y=[]
        #for i in lines:
            #if i[0] in chars:
             #   continue
            #else:
                #x.append(i.split()[0])
                #y.append(i.split()[0])
        x=[float(x.split()[0]) for x in lines if x[0] not in chars ]
        y=[float(x.split()[1]) for x in lines if x[0] not in chars ]
        
        x_corr=[j+5.14 for j in x]
        y_scaled=[i-min(y) for i in y]
        #frcs_tmp=[[float(x.split()[i]) for x in lines[1:]] for i in range(6)]
        data.append([x_corr,y_scaled])
       
        names.append(file.replace("/Volumes/My Passport/RATIO/RATIO/WP4/FFs/umbrella/MOD-FRC/BIG/","").replace(
            "/umbrella_30_200/profile.xvg",''
        ))
for file in histo:
     with open(file,'r') as ifile:
        lines=ifile.readlines()
        hist_x=[float(x.split()[0]) for x in lines if x[0] not in chars]
        hist_y=[[float(x.split()[i]) for x in lines if x[0] not in chars] for i in range(1,26)]
        hist_x_corr=[j+5.14 for j in hist_x]
        data_hist.append([hist_x_corr,hist_y])
    
for file in dens:
    with open(file,'r') as ifile:
        lines=ifile.readlines()
        x=[float(x.split()[0]) for x in lines if x[0] not in chars ]
        y=[float(x.split()[1]) for x in lines if x[0] not in chars ]
        rho.append([x,y])


#%%
labelsX="z / nm"# / eV","frcY_QM / eV","frcZ_QM / eV"]
labelsY="E / $KJ\ mol^{-1}$"# / eV","frcY_ML / eV","frcZ_ML / eV"]
limx=[6,16]
limy=[-1,40]
fig = plt.figure(figsize=(20,10)) 
#print(ftot[1])
#rint(frcs[0])
for k in range(6):
        ax=fig.add_subplot(2,3,k+1)
        ax.plot(data[k][0],data[k][1])
        ax.axvline(x=10.6, color = '#606060')
        ax.axvline(x=8.4, color = '#606060')        
        ax.set_xlim(limx)
        ax.set_ylim(limy)
        plt.xlabel(labelsX)
        plt.ylabel(labelsY)
        plt.title(labs[names[k]])
        #ax.text(.99,.02,"AvgE= {:.2e} eV".format(avg_err[k][i]), ha='right',transform=ax.transAxes)

    #plt.setp()
plt.rcParams.update({'font.size': 18})
plt.tight_layout()
plt.show()
#%%
labelsX_hist="z / nm"# / eV","frcY_QM / eV","frcZ_QM / eV"]
labelsY_hist="count"# / eV","frcY_ML / eV","frcZ_ML / eV"]
limx=[6.5,15]
limy=[-1,350000]
fig = plt.figure(figsize=(20,10)) 
#print(ftot[1])
#rint(frcs[0])
#print(data_hist[0][2][1])
#print(data_hist[5][1][0])
for k in range(6):
        ax=fig.add_subplot(2,3,k+1)
        for i in range(0,len(data_hist[k][1])):
            ax.plot(data_hist[k][0],data_hist[k][1][i])
        ax.set_xlim(limx)
        #ax.set_ylim(limy)
        plt.xlabel(labelsX_hist)
        plt.ylabel(labelsY_hist)
        plt.title(labs[names[k]])
        #ax.text(.99,.02,"AvgE= {:.2e} eV".format(avg_err[k][i]), ha='right',transform=ax.transAxes)
plt.rcParams.update({'font.size': 18})
    #plt.setp()
plt.tight_layout()
plt.show()

#%%
fig2 = plt.figure(figsize=(10,5),dpi=150) 
clrs=["r","b"]
plt.axvline(x=10.6, color = '#606060')
plt.axvline(x=8.4, color = '#606060')
for i,item in enumerate(rho):
    plt.plot(item[0],item[1],color=clrs[i])


    
plt.xlabel(r"$z\ /\ nm$")
plt.ylabel(r'$\rho\ /\ kg\ m^{-3}$')
plt.legend(["DES","THF"])
plt.show()

# %%
