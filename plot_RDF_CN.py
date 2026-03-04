import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Extract data from file')
parser.add_argument('--data', type=str, help='Input data file', required=True)
parser.add_argument('--output', type=str, help='Output file', required=True)
parser.add_argument('--MM', type=int, help='Denominator exponent', required=True)
parser.add_argument('--NN', type=int, help='Numerator exponent', required=True)
parser.add_argument('--r0', type=float, help='r0 value', default=1.0)
parser.add_argument('--color', type=int, help='Plot color', default=0)
parser.add_argument('--title', type=str, help='Plot title', default="")
args= parser.parse_args()


font = {'family' : 'Formular',
        'weight' : 'normal',
        'size'   : 46}
mpl.rc('font', **font)
mpl.rcParams['axes.linewidth'] = 3
mpl.rcParams['lines.linewidth'] = 3


pres_colors=["#c1272d","#0000a7","#eba938","#008176","#b3b3b3","#4cb944"]

data=np.loadtxt(args.data, unpack=True,comments=['@','&'])

CN_X=data[0]
CN_Y=[(1-(r/ args.r0)**args.NN)/(1-(r/args.r0)**args.MM) for r in data[0]]
#CN_Y=[x if x is not nan else 0.5 for x in CN_Y]
for i in CN_Y:
    if np.isnan(i):
        CN_Y[CN_Y.index(i)]=0.5

fig=plt.figure(figsize=(16,12))
label=f"RDF {args.data.split('.')[0].split('_')[0]}-{args.data.split('.')[0].split('_')[1]}"
plt.plot(data[0], data[1]/np.max(data[1]), label=label,color=pres_colors[args.color])
plt.plot(CN_X, CN_Y, label='CN Function',color='black')
plt.scatter([args.r0], [0.5], marker='D', color='black',s=50 )
plt.vlines(args.r0,0,0.5,colors='black',linestyles='dashed',label=f'$r_0$ = {args.r0} Å')

plt.xlabel('Distance (Angstrom)')
plt.ylabel('Normalized RDF')
plt.title(args.title)
plt.legend()
plt.tight_layout()
plt.savefig(args.output, format='png')
plt.close()

