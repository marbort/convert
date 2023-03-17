import matplotlib.pyplot as plt
import sys, os
import glob
import matplotlib.colors as mcolors

files=["test.f.out","test2.f.out"]
frcs=[]
names=[]
for file in files:
    with open(file,'r') as ifile:
        lines=ifile.readlines()
        frcs_tmp=[[float(x.split()[i]) for x in lines[1:]] for i in range(6)]
        frcs.append(frcs_tmp)
        names.append(file)

print(len(frcs[0][0]))
err=[[[abs(frcs[i][j][k]-frcs[i][j+3][k]) for k in range(len(frcs[i][j])) ]for j in range(3)  ] for i in range(len(frcs))]
avg_err=[[sum(x)/len(x) for x in l] for l in err]
print("{:.2e}".format(avg_err[0][2]))
print(mcolors.CSS4_COLORS[list(mcolors.CSS4_COLORS.keys())[0]])


labelsX=["frcX_QM / eV","frcY_QM / eV","frcZ_QM / eV"]
labelsY=["frcX_ML / eV","frcY_ML / eV","frcZ_ML / eV"]
fig = plt.figure() 
#print(ftot[1])
#rint(frcs[0])
for k in range(len(frcs)):
    for i in range(3):
        ax=fig.add_subplot(len(frcs),3,3*k+i+1)
        ax.plot(frcs[k][0+i],frcs[k][3+i],'x',markersize=5)
        plt.xlabel(labelsX[i])
        plt.ylabel(labelsY[i])
        plt.title(names[k])
        ax.text(.99,.02,"AvgE= {:.2e} eV".format(avg_err[k][i]), ha='right',transform=ax.transAxes)
                 
    #plt.setp()
plt.tight_layout()
plt.show()

