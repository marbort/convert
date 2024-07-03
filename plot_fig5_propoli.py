import pandas as pd
import numpy
import matplotlib.pyplot as plt
import sys
import matplotlib as mpl
import numpy as np

font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : 26}
mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['lines.linewidth'] = 3


classes={
    "chalchones and dihydrochalcones": [1,5],
    "flavanones": [6,14],
    "flavones"  : [15,18],
    "caffeic acid derivatives" : [19,28] 
}

classes_labels=["chalchones and\ndihydrochalcones","flavanones","flavones","caffeic acid\nderivatives"]

fractions=["f4","f5","f6"]

averages={"HO":["HO"],"CH3O":["CH3O"],"ROO":["HOO","CH3OO","CH2CHOO","CH2CHCHCHOO"],"L":["L31","L32","L33","L34"]}
average_data={"HO":{},"CH3O":{},"ROO":{},"L":{}}
average_data_classes={}
rad_title=["HO^•","CH_3O^•","ROO^•","L^•"]

blues = mpl.cm.get_cmap('Blues')
greens = mpl.cm.get_cmap('Greens')
reds = mpl.cm.get_cmap('Reds')
yellows = mpl.cm.get_cmap('Wistia')
colors_chal = [greens(x) for x in np.linspace(0,1,5)]
colors_flava = [reds(x) for x in np.linspace(0,1,9)]
colors_flavo = [yellows(x) for x in np.linspace(0.2,0.8,4)]
colors_caffe = [blues(x) for x in np.linspace(0,1,10)]
colors=colors_chal+colors_flava+colors_flavo+colors_caffe

colors_cls=['yellowgreen','orangered','gold','navy']
print(colors)
    
for j in averages.keys():
        average_data_classes[j]={}
        for i in classes.keys():
            average_data_classes[j][i]={}
    
data=pd.read_excel(sys.argv[1])

for index,row in data.iterrows():
    for avg in averages.keys():
        average_data[avg][int(row['ros'])]={}
        for fraction in fractions:
            average_data[avg][int(row['ros'])][fraction]=sum([row[x] for x in averages[avg]])*row[f"{fraction} %tot"]/len(averages[avg])/100
for rad in average_data_classes:
    for cls in average_data_classes[rad]:
        for fraction in fractions:
            average_data_classes[rad][cls][fraction]=sum([average_data[rad][x][fraction] for x in range(classes[cls][0],classes[cls][1]+1)])
#print(average_data)            
#print(average_data_classes)

fig=plt.figure(figsize=(15,20))

for i,rad in enumerate(averages):
    plt.subplot(len(averages),2,2*i+1)
    bottom=np.zeros(len(fractions))
    for j,cls in enumerate(average_data_classes[rad]):
        values=[average_data_classes[rad][cls][x] for x in average_data_classes[rad][cls]]
        plt.bar([x.upper() for x in fractions],values,color=colors_cls[j],bottom=bottom)
        bottom+=values
    plt.title("$\mathrm{\mathsf{"+f"{rad_title[i]}"+"}}$")
    
    plt.ylabel("$\mathrm{\mathsf{\Delta G^{0\prime}_{HAT}\ (kcal/mol)}}$")
    plt.ylim(-20,4)
    if i==len(averages)-1:
        plt.xlabel("Fraction")
        plt.legend(classes_labels,loc=(0.02,0.02),fontsize=16,ncols=2)
    
        
    
    plt.subplot(len(averages),2,2*i+2)
    bottom2=np.zeros(len(fractions))
    for j,species in enumerate(average_data[rad]):
            values2=[average_data[rad][species][x] for x in average_data[rad][species]]
            plt.bar([x.upper() for x in fractions],values2,color=colors[j],bottom=bottom2)
            bottom2+=values2
    plt.title("$\mathrm{\mathsf{"+f"{rad_title[i]}"+"}}$")
    
    plt.ylim(-19,4)
    if i==len(averages)-1:
        plt.xlabel("Fraction")
        plt.legend(average_data[rad].keys(),loc=(0.02,0.02),ncols=5,fontsize=14)





    
plt.tight_layout()
#plt.show()
plt.savefig('Fig5.png',format="png",dpi=300)
        





