import numpy as np
import argparse 
import matplotlib.pyplot as plt
import os
import matplotlib as mpl

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    print(array[idx])
    return array[idx]

def extract_data(input,split,min,fac):
    data={}
    freesplit=[]
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    data['full']={'cv':[float(x.split()[0])*fac for x in lines if "#" not in x],'free':[float(x.split()[1]) for x in lines if "#" not in x]}
    with open("{}_equil".format(input),'r') as ifile:
        lines=ifile.readlines()
    data['equil']={'cv':[float(x.split()[0])*fac for x in lines if "#" not in x],'free':[float(x.split()[1]) for x in lines if "#" not in x]}
    for i in range(split):
        with open("{}_{}".format(input,i),'r') as ifile:
            lines=ifile.readlines()
        data["split_{}".format(i)]={"cv":[],"free":[]}
        for line in lines:
            if '#' not in line:
                #if 'inf' not in line:
                    data["split_{}".format(i)]['cv'].append(float(line.split()[0])*fac)
                    data["split_{}".format(i)]['free'].append(float(line.split()[1]))
        if i != 0:
            freesplit.append(data["split_{}".format(i)]['free'])                 
    data['mean_free']=[np.mean(np.ma.masked_invalid(x)) for x in zip(*freesplit)]
    if min:
        cv_min=find_nearest(data['split_0']['cv'],min)
        data['mean_free_scaled']=[x-data['mean_free'][data['split_0']['cv'].index(cv_min)] for x in data['mean_free']]
        data['full']['free_scaled']=[x-data['full']['free'][data['split_0']['cv'].index(cv_min)] for x in data['full']['free']]
        data['equil']['free_scaled']=[x-data['equil']['free'][data['split_0']['cv'].index(cv_min)] for x in data['equil']['free']]
    data['std_free']=[np.std(np.ma.masked_invalid(x)) for x in zip(*freesplit)]
    data['std_free_mean']=[np.std(np.ma.masked_invalid(x))/np.sqrt(split-1) for x in zip(*freesplit)]
    #print(data['mean_free_scaled'][data['split_0']['cv'].index(cv_min)])
    print(data['std_free'][12],data['std_free_mean'][12])
    
    return(data)

def plot_data(data,min,full,equil,input,kcal,x,y,labels,inverted,pres,wd=9,hg=9):
    fig=plt.figure(figsize=(wd,hg),dpi=150)
    font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : 16}
    mpl.rc('font', **font)
    """
    for j,i in enumerate(data):
        plt.errorbar([x*fac[j] for x in data[i]['cv']],data[i]['free'],yerr=data[i]['err'],ecolor=erclr,marker='o',markersize='3', label=i)
    plt.xlabel('CV')
    plt.ylabel('Free Energy ($kJ\ mol^{-1}$)')
    plt.legend()
    plt.savefig('umbrella_fes.png',format='png')
    """
    kj_to_kcal=0.239006
    erclr=["1f77b4", "ff7f0e", "2ca02c", 'd62728', '9467bd', '8c564b', 'e377c2', '7f7f7f', 'bcbd22', '17becf']
    erclr_rgba=[[int("".join(x[0:2]),16)/255,int("".join(x[2:4]),16)/255,int("".join(x[4:6]),16)/255,0.2] for x in erclr]
    pres_colors=["2b38ff","f7059b","17d9ff","000000","4cb944"]
    pres_colors_rgba=[[int("".join(x[0:2]),16)/255,int("".join(x[2:4]),16)/255,int("".join(x[4:6]),16)/255,1] for x in pres_colors]
    for j,i in enumerate(data):
        if pres:
            clr=pres_colors_rgba[j]
        else:
            clr="C{}".format(j)
        if equil:
            if kcal:
                if min:
                   
                    markers,caps,bars=plt.errorbar([x for x in data[i]["equil"]['cv']],[x*kj_to_kcal for x in data[i]['equil']['free_scaled']],yerr=[x*1.96*kj_to_kcal if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']],
                            ecolor=erclr_rgba[j], marker='o',markersize='3', label="_".join([i.split('/')[x] for x in range(-6,-1)]),color=clr) #ecolor=erclr_rgba
                    plt.fill_between(data[i]['equil']['cv'],[x*kj_to_kcal-[x*1.96*kj_to_kcal if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']][k] for k,x in enumerate(data[i]['equil']['free_scaled'])],
                                     [x*kj_to_kcal+[x*1.96*kj_to_kcal if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']][k] for k,x in enumerate(data[i]['equil']['free_scaled'])],linewidth=0,alpha=0.1,color=clr)
                    if full:
                        markers,caps,bars=plt.plot(data[i]['full']['cv'],[x*kj_to_kcal for x in data[i]['full']['free_scaled']],'--',color=clr)
                else:
                    markers,caps,bars=plt.errorbar([x for x in data[i]["equil"]['cv']],[x*kj_to_kcal for x in data[i]['equil']['free']],yerr=[x*1.96*kj_to_kcal if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']],
                                ecolor=erclr_rgba[j],marker='o',markersize='3', label="_".join([i.split('/')[x] for x in range(-6,-1)]),color=clr) #ecolor=erclr_rgba,
                    plt.fill_between(data[i]['equil']['cv'],[x*kj_to_kcal-[x*1.96*kj_to_kcal if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']][k] for k,x in enumerate(data[i]['equil']['free'])],
                                     [x*kj_to_kcal+[x*1.96*kj_to_kcal if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']][k] for k,x in enumerate(data[i]['equil']['free'])],linewidth=0,alpha=0.1,color=clr)
                    if full:
                        markers,caps,bars=plt.plot(data[i]['full']['cv'],[x*kj_to_kcal for x in data[i]['full']['free']],'--',color=clr)
            else:
                if min:
                    #print(data[i]['equil']['free_scaled']) 
                    markers,caps,bars=plt.errorbar([x for x in data[i]["equil"]['cv']],data[i]['equil']['free_scaled'],yerr=[x*1.96 if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']],
                            ecolor=erclr_rgba[j], marker='o',markersize='3', label="_".join([i.split('/')[x] for x in range(-6,-1)]),color=clr) #ecolor=erclr_rgba
                    plt.fill_between(data[i]['equil']['cv'],[x-[x*1.96 if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']][k] for k,x in enumerate(data[i]['equil']['free_scaled'])],
                                     [x+[x*1.96 if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']][k] for k,x in enumerate(data[i]['equil']['free_scaled'])],linewidth=0,alpha=0.1,color=clr)
                    if full:
                        markers,caps,bars=plt.plot(data[i]['full']['cv'],data[i]['full']['free_scaled'],'--',color=clr)
                        
                else:
                    markers,caps,bars=plt.errorbar([x for x in data[i]["equil"]['cv']],data[i]['equil']['free'],yerr=[x*1.96 if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']],
                                ecolor=erclr_rgba[j],marker='o',markersize='3', label="_".join([i.split('/')[x] for x in range(-6,-1)]),color=clr) #ecolor=erclr_rgba,
                    plt.fill_between(data[i]['equil']['cv'],[x-[x*1.96 if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']][k] for k,x in enumerate(data[i]['equil']['free'])],
                                     [x+[x*1.96 if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']][k] for k,x in enumerate(data[i]['equil']['free'])],linewidth=0,alpha=0.1,color=clr)
                    if full:
                        markers,caps,bars=plt.plot(data[i]['full']['cv'],data[i]['full']['free'],'--',color=clr)
        else:
            if kcal:
                if min:
                    markers,caps,bars=plt.errorbar([x for x in data[i]["split_0"]['cv']],[x*kj_to_kcal for x in data[i]['mean_free_scaled']],yerr=[x*1.96*kj_to_kcal if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']],
                            ecolor=erclr_rgba[j], marker='o',markersize='3', label="_".join([i.split('/')[x] for x in range(-6,-1)]),color=clr) #ecolor=erclr_rgba
                    if full:
                        plt.plot(data[i]['full']['cv'],[x*kj_to_kcal for x in data[i]['full']['free_scaled']],'--',color=clr)
                else:
                    markers,caps,bars=plt.errorbar([x for x in data[i]["split_0"]['cv']],[x*kj_to_kcal for x in data[i]['mean_free']],yerr=[x*1.96*kj_to_kcal if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']],
                                ecolor=erclr_rgba[j],marker='o',markersize='3', label="_".join([i.split('/')[x] for x in range(-6,-1)]),color=clr) #ecolor=erclr_rgba,
                    if full:
                        markers,caps,bars=plt.plot(data[i]['full']['cv'],[x*kj_to_kcal for x in data[i]['full']['free']],'--',color=clr)
            else:
                if min:
                    
                    markers,caps,bars=plt.errorbar([x for x in data[i]["split_0"]['cv']],data[i]['mean_free_scaled'],yerr=[x*1.96 if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']],
                            ecolor=erclr_rgba[j], marker='o',markersize='3', label="_".join([i.split('/')[x] for x in range(-6,-1)]),color=clr) #ecolor=erclr_rgba
                    if full:
                        markers,caps,bars=plt.plot(data[i]['full']['cv'],data[i]['full']['free_scaled'],'--',color=clr)
                else:
                    markers,caps,bars=plt.errorbar([x for x in data[i]["split_0"]['cv']],data[i]['mean_free'],yerr=[x*1.96 if x is not np.ma.masked else 0 for x in data[i]['std_free_mean']],
                                ecolor=erclr_rgba[j],marker='o',markersize='3', label="_".join([i.split('/')[x] for x in range(-6,-1)]),color=clr) #ecolor=erclr_rgba,
                    if full:
                        markers,caps,bars=plt.plot(data[i]['full']['cv'],data[i]['full']['free'],'--',color=clr)
        [bar.set_alpha(0.0) for bar in bars]
    #plt.plot([x for x in data[i]["split_0"]['cv']],np.zeros(len(data[i]["split_0"]['cv'])))
    plt.xlabel(x)
    plt.ylabel(y)
    #plt.xlim([0.8,2.1])
    #plt.ylim([-5,26])
    
    if labels:
        #leg=plt.legend(labels,loc='upper center')
        leg=plt.legend(labels)
        for lh in leg.legendHandles: 
            lh.set_alpha(1)
    #else:
    #    plt.legend()
    if inverted:
        plt.gca().invert_xaxis()
    plt.savefig('umbrella_fes_{}.png'.format(os.path.basename(input[0])),format='png')

    

def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--input', dest='input', 
                        type=str, help='input data',nargs="+")
    parser.add_argument('--fac', dest='fac', 
                        type=float, help='cv conversion factors',nargs="+")
    parser.add_argument('--split', dest='split', 
                        type=int, help='number of split for each umbrella',nargs="+")
    parser.add_argument('--min', dest='min', 
                        type=float, help='CV value for 0')
    parser.add_argument('--full', dest='full', action='store_true',
                         help='plot also FES using whole trajectory')
    parser.add_argument('--equil', dest='equil', action='store_true',
                         help='plot also FES using all but first block')
    parser.add_argument('--kcal', dest='kcal', action='store_true',
                         help='plot FES in kcal')
    parser.add_argument('--xlabel', dest='xlabel', default='CV', type=str,
                         help='x label')
    parser.add_argument('--ylabel', dest='ylabel', default='Free Energy ($kJ\ mol^{-1}$)', type=str,
                         help='y label')
    parser.add_argument('--labels', dest='labels', default=None, type=str, nargs='+',
                         help='Plot labels label')
    parser.add_argument('--inverted', dest='inverted', action='store_true',
                         help='Inverts X axis')
    parser.add_argument('--pres', dest='pres', action='store_true',
                         help='use presentation colors')
    args = parser.parse_args()
    data={}
    for i,input in enumerate(args.input):
        data[input]=extract_data(input,args.split[i],args.min,args.fac[i])
    #print(data)
    plot_data(data,args.min,args.full,args.equil,args.input,args.kcal,args.xlabel,args.ylabel,args.labels,args.inverted,args.pres)
    
    #print(len(data['free']),len(data['err']))
    
if __name__=="__main__":
    main()