#!/Users/mbortoli/opt/anaconda3/bin/python

import sys, os
import regex as re
from atom_dict import Atomic_number
import argparse
import cclib
import numpy as np
import glob
import subprocess

#arg1 filename; arg2 CI/noCI; arg3 numstruct; arg4 SI/XMOL


def split_g16_calc(input):
    splits=[]
    tmpfiles=[]
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    split_idx=[0]+[i for i,line in enumerate(lines) if "Normal" in line]
    splits=[lines[split_idx[x-1]+1:split_idx[x]+1] for x in range(1,len(split_idx))]
    opt_freq=[[1 for line in split if "Guess=TCheck" in line] for split in splits]
    opt_freq_num=[i for i,num in enumerate(opt_freq) if len(num)>0]
    splits_final=[splits[i-1]+splits[i] for i in opt_freq_num]
    for i,split in enumerate(splits):
        if i in opt_freq_num or i in [x-1 for x in opt_freq_num]:
            pass
        else:
            splits_final.append(splits[i])
    solvents=[[x.split()[2][:-1] for x in split if "Solvent              :" in x] for split in splits ]
    solvents_final=[solvents[i-1]+solvents[i] for i in opt_freq_num]
    for i,solvent in enumerate(solvents):
        if i in opt_freq_num or i in [x-1 for x in opt_freq_num]:
            pass
        else:
            solvents_final.append(solvents[i])
    solv_models=[[x.split()[3][:-1] for x in split if "Atomic radii         :" in x] for split in splits ]
    solv_models_final=[solv_models[i-1]+solv_models[i] for i in opt_freq_num]
    for i,solv_model in enumerate(solv_models):
        if i in opt_freq_num or i in [x-1 for x in opt_freq_num]:
            pass
        else:
            solv_models_final.append(solv_models[i])
    for i,split in enumerate(splits_final):
        with open(input+'tmp'+str(i),'w') as ofile:
            tmpfiles.append(input+'tmp'+str(i))
            for line in split:
                ofile.write(line)
    return(tmpfiles,solvents_final,solv_models_final)
   
        
def parse_files(input,solvents,solvent_models,coords,tipo):
    for j, file in enumerate(input):
        #if file.endswith(args.ext):
        try:
            with open(file,'r') as ifile:
                    data=cclib.io.ccread(ifile)
                    basname=re.sub(r'tmp[0-9]','',os.path.splitext(ifile.name)[0])
                    elm=[]
                    for i in data.atomnos:
                        elm.append(Atomic_number[i])
                    with open(basname+'.xyz','w') as ofile:
                        fullcrd=list(zip(elm,data.atomcoords[-1]))
                        #print(fullcrd[0][1][2])
                        ofile.write('{}\n'.format(data.natom))
                        ofile.write('charge={}, mult={}, Total Energy={:.5f} Ha\n'.format(data.charge,data.mult,data.scfenergies[-1]/27.21138505))   
                        for i in range(len(fullcrd)):
                            ofile.write('{}\t'.format(fullcrd[i][0]))
                            for k in range(len(fullcrd[i][1])):
                                ofile.write('{:8.5f}\t'.format(fullcrd[i][1][k]))
                            ofile.write('\n')
                    if tipo == "SI":
                        with open("SI.txt",'a') as sifile:
                            sifile.write(basname+'\n')
                            sifile.write(f'{data.natom}\n')
                            sifile.write(f"Basis Set {data.metadata['basis_set']}\n")
                            sifile.write(f"Multiplicity {data.mult}\n")
                            sifile.write(f'Solvent Radii {solvent_models[j]}\n')
                            sifile.write(f'Solvent {solvents[j]}\n')
                            sifile.write(f"Electronic Energy={data.scfenergies[-1]/27.21138505:.5f} Ha\n")
                            try:
                                sifile.write(f"Sum of electronic and thermal enthalpies: {data.enthalpy:5f} Ha\n")
                                sifile.write(f"Sum of electronic and thermal free energy: {data.freeenergy:5f} Ha\n")
                                sifile.write(f"Lowest Frequency: {data.vibfreqs[0]:3f} 1/cm\n")
                            except:
                                pass
                            if coords[j]==1:    
                                for i in range(len(fullcrd)):
                                    sifile.write(f'{fullcrd[i][0]:4s}  ')
                                    for k in range(len(fullcrd[i][1])):
                                        sifile.write('{:8.5f}    '.format(fullcrd[i][1][k]))
                                    sifile.write('\n')
                                sifile.write('\n')
                            else:
                                sifile.write("\nSee previous structure\n\n")
                #np.savetxt(ofile,crds)
            os.remove(file)
            #print('ciao')
        except:
                print("Data not found skipping file {}".format(file))        


def main():

    parser = argparse.ArgumentParser(description='Exctract xyz coordinates')
    #parser.add_argument('--CI', dest='CI', default=False, help='CI output not converged present in folder')
    parser.add_argument('--type', dest='tipo', default='XMOL',  help='Type of print: SI for complete SI; XMOL for separate xyz files')
    parser.add_argument('--ext', dest='ext', default='.log', help='Extension of the Gaussian results file. Default .log')
    parser.add_argument('--recursive', dest='recursive', default=False, action='store_true', help='Look recursively in all folders in current dir')
    parser.add_argument('--coords',dest='coords',type=int,nargs='+',help='Array to print coordinates of stage if 1. If 0 do not print. Must be the same \
                        length of number of calculations (opt+freq counts as one)')
    args = parser.parse_args()

    if args.recursive:
        files=glob.glob('./**/*'+args.ext,recursive=True)
    else:
        files=sorted(glob.glob("./*"+args.ext))
    
    
    print(files)
    for file in files:
        tmpfiles,solvents,solv_models=split_g16_calc(file)
        #print(tmpfiles,solvents)
        parse_files(tmpfiles,solvents,solv_models,args.coords,args.tipo)
    
   
    


if __name__=="__main__":
    main()

