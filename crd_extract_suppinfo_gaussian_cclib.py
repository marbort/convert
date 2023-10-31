#!/Users/mbortoli/opt/anaconda3/bin/python

import sys, os
import regex as re
from atom_dict import Atomic_number
import argparse
import cclib
import numpy as np

#arg1 filename; arg2 CI/noCI; arg3 numstruct; arg4 SI/XMOL


parser = argparse.ArgumentParser(description='Exctract xyz coordinates')
#parser.add_argument('--CI', dest='CI', default=False, help='CI output not converged present in folder')
parser.add_argument('--type', dest='tipo', default='XMOL',  help='Type of print: SI for complete SI; XMOL for separate xyz files')
parser.add_argument('--ext', dest='ext', default='.log', help='Extension of the Gaussian results file. Default .log')
args = parser.parse_args()


for file in sorted(os.listdir('./')):
    if file.endswith(args.ext):
        try:
           with open(file,'r') as ifile:
               data=cclib.io.ccread(ifile)
               basname=os.path.splitext(ifile.name)[0]
               elm=[]
               for i in data.atomnos:
                   elm.append(Atomic_number[i])
               with open(basname+'.xyz','w') as ofile:
                   fullcrd=list(zip(elm,data.atomcoords[-1]))
                   #print(fullcrd[0][1][2])
                   ofile.write('{}\n'.format(data.natom))
                   ofile.write('charge={}, mult={}, Total Energy={:.5f} Ha\n'.format(data.charge,data.mult,data.scfenergies[-1]/27.2114))   
                   for i in range(len(fullcrd)):
                       ofile.write('{}\t'.format(fullcrd[i][0]))
                       for k in range(len(fullcrd[i][1])):
                           ofile.write('{:8.5f}\t'.format(fullcrd[i][1][k]))
                       ofile.write('\n')
               if args.tipo == "SI":
                   with open("SI.txt",'a') as sifile:
                       sifile.write(basname+'\n')
                       #sifile.write('{}\n'.format(data.natom))
                       sifile.write('Total Energy={:.5f} Ha\nNimag=0 \n'.format(data.scfenergies[-1]/27.2114))   
                       for i in range(len(fullcrd)):
                           sifile.write(f'{fullcrd[i][0]:4s}  ')
                           for k in range(len(fullcrd[i][1])):
                               sifile.write('{:8.5f}    '.format(fullcrd[i][1][k]))
                           sifile.write('\n')
                       sifile.write('\n')
                #np.savetxt(ofile,crds)
        except:
              print("Data not found skipping file {}".format(file))        


