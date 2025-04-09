
import argparse
import glob
import os
import subprocess
import sys
import numpy as np

def skip_lines(file,lines,skip):
    with open(f"{file.name}_skipped",'w') as ofile:
        for line in lines[skip:]:
            ofile.write(line)

def calc_min_max(file_min,file_max):
    data_min=np.loadtxt(file_min,unpack=True)
    data_max=np.loadtxt(file_max,unpack=True)
    return(np.min(data_min[1]),np.max(data_max[1]))
    

def extract_params(file,pos,sep,restart,fac):
    at=file.split('-')[0].split(f'{sep}')[pos]
    if restart:
       cp2k_input=f"cp2k_{at}_rst.inp"
       
    else:
        cp2k_input=f"cp2k_{at}.inp"
    try:
        with open(cp2k_input,'r') as cpinp:
            lines=cpinp.readlines()
            for line in lines:
                if 'K [' in line:
                    kappa=line.split()[-1]
        
    except:
        print(f"Skipping {at} restart. File not found")
        kappa=0
    return(float(kappa)*2,round(float(at)/fac,2))

def extract_params_restart(file,pos,sep,idx,fac):
    at=file.split('-')[0].split(f'{sep}')[pos]
    cp2k_restart=file.replace(f"-COLVAR.metadynLog_skipped_{idx}",'-1.restart')
    kappa=0
    try:
        with open(cp2k_restart,'r') as cpinp:
            lines=cpinp.readlines()
        for line in lines:
            if 'K  ' in line:
                kappa=line.split()[-1]
        return(round(float(kappa)*2*2625.5002,1),round(float(at)/fac,2)) 
    except Exception as X:
        print(X)
        print(f"Skipping {at} restart. File not found")
        return(round(float(kappa)*2*2625.5002,1),round(float(at)/fac,2))
    
    

def split_file(file,lines,parts):
    splits=len(lines)//parts
    for i in range(parts):
        with open(f"{file}_{i}",'w') as ofile:
            for line in lines[splits*i:splits*(i+1)]:
                ofile.write(line)
            if i == parts-1:
                try:
                    ofile.write(lines[-1])
                except:
                    os.remove(ofile.name)
                   
def do_wham(parts,start,end,nbins,temp,tol,periodic):
    if periodic is not None:
        if periodic=="pi":    
            subprocess.run(["/home/marco/WHAM/wham/wham/wham", "Ppi", start, end, nbins, tol, temp, "0", "sims.txt", "out.dat_all"])
            subprocess.run(["/home/marco/WHAM/wham/wham/wham", "Ppi", start, end, nbins, tol, temp, "0", "sims_equil.txt", "out.dat_all_equil"])
            subprocess.run(["/home/marco/WHAM/wham/wham/wham", "Ppi", start, end, nbins, tol, temp, "0", "sims_equil.txt", "out.dat_all_equil_errors","25", f"{np.random.randint(0,1000)}"])
            for i in range(parts):
                subprocess.run(["/home/marco/WHAM/wham/wham/wham", "Ppi", start, end, nbins, tol, temp, "0", f"sims_{i}.txt", f"out.dat_all_{i}"])
        elif periodic=="360":    
            subprocess.run(["/home/marco/WHAM/wham/wham/wham", "P", start, end, nbins, tol, temp, "0", "sims.txt", "out.dat_all"])
            subprocess.run(["/home/marco/WHAM/wham/wham/wham", "P", start, end, nbins, tol, temp, "0", "sims_equil.txt", "out.dat_all_equil"])
            subprocess.run(["/home/marco/WHAM/wham/wham/wham", "P", start, end, nbins, tol, temp, "0", "sims_equil.txt", "out.dat_all_equil_errors", "25", f"{np.random.randint(0,1000)}"])
            for i in range(parts):
                subprocess.run(["/home/marco/WHAM/wham/wham/wham", "Ppi", start, end, nbins, tol, temp, "0", f"sims_{i}.txt", f"out.dat_all_{i}"])
        else:
            subprocess.run(["/home/marco/WHAM/wham/wham/wham", f"P{periodic}", start, end, nbins, tol, temp, "0", "sims.txt", "out.dat_all"])
            subprocess.run(["/home/marco/WHAM/wham/wham/wham", f"P{periodic}", start, end, nbins, tol, temp, "0", "sims_equil.txt", "out.dat_all_equil"])
            subprocess.run(["/home/marco/WHAM/wham/wham/wham", f"P{periodic}", start, end, nbins, tol, temp, "0", "sims_equil.txt", "out.dat_all_equil_errors", "25", f"{np.random.randint(0,1000)}"])
            for i in range(parts):
                subprocess.run(["/home/marco/WHAM/wham/wham/wham", f"P{periodic}", start, end, nbins, tol, temp, "0", f"sims_{i}.txt", f"out.dat_all_{i}"])
    else:
        subprocess.run(["/home/marco/WHAM/wham/wham/wham", start, end, nbins, tol, temp, "0", "sims.txt", "out.dat_all"])
        subprocess.run(["/home/marco/WHAM/wham/wham/wham", start, end, nbins, tol, temp, "0", "sims_equil.txt", "out.dat_all_equil"])
        subprocess.run(["/home/marco/WHAM/wham/wham/wham", start, end, nbins, tol, temp, "0", "sims_equil.txt", "out.dat_all_equil_errors", "25", f"{np.random.randint(0,1000)}"])
        for i in range(parts):
            subprocess.run(["/home/marco/WHAM/wham/wham/wham", start, end, nbins, tol, temp, "0", f"sims_{i}.txt", f"out.dat_all_{i}"])

    
    


def main():
    
    parser = argparse.ArgumentParser(description='Plot data')

    parser.add_argument('--parts' , dest='parts',default=5, type=int,help='Number of splits')
    parser.add_argument('--pos' , dest='pos',default=4,type=int,help='Restrain position in filename')
    parser.add_argument('--sep' , dest='sep',default='_',help='Separator to get restrain position')
    parser.add_argument('--restart' , dest='restart',action='store_true', default=False, help='Use CP2K restart for restraint parameters',)
    parser.add_argument('--skip' , dest='skip',default=1000, type=int,help='Number of lines to skip from colvar file')
    parser.add_argument('--nbins' , dest='nbins', type=str, default="0.1",help='Size of the bins for WHAM')
    parser.add_argument('--temp' , dest='temp',type=str, default="300",help='Temperature for WHAM')
    parser.add_argument('--fac' , dest='fac',type=float, default=1.89,help='Conversion factor between at in filename and in colvar. Default: 1.89')
    parser.add_argument('--tol' , dest='tol',type=str, default="1e-5",help='Tolerance for WHAM. Default 1e-5')
    parser.add_argument('--periodic' , dest='periodic',type=str, default=None,help='Periodicity of the CV. Default None')
    parser.add_argument('--wham' , dest='wham',action='store_true', default=False, help='Use WHAM range')
    
    args = parser.parse_args()
    
    print(f"Conversion factor for CV between files and colvar: {args.fac}")
    
    files_all=glob.glob('*metadynLog')
    files=[x for x in files_all if "equil" not in x]
    files_extremes=[sorted(files)[0],sorted(files)[-1]]
    minwham,maxwham=calc_min_max(files_extremes[0],files_extremes[1])
    start=files_extremes[0].split('-')[0].split(f'{args.sep}')[args.pos]
    end=files_extremes[1].split('-')[0].split(f'{args.sep}')[args.pos]
    print(f"Starting from {start} to {end}")
    print(f"WHAM full range: {minwham*args.fac:.2f} to {maxwham*args.fac:.2f}")
    print(f"WHAM allowed range: {float(start)} to {float(end)}")
    if not args.wham:
        if minwham*args.fac < float(start):
            minwham=float(start)/args.fac   
        if maxwham*args.fac > float(end):
            maxwham=float(end)/args.fac
    print(f"Get range from trajectory: {args.wham}")
    print(f"WHAM range: {minwham*args.fac:.2f} to {maxwham*args.fac:.2f}")
    nbins=int((maxwham-minwham)*args.fac/float(args.nbins))
    
    for file in files:
        with open(file,'r') as ifile:
            lines=ifile.readlines()
            skip_lines(ifile,lines,args.skip)
            split_file(file,lines,args.parts)
    files_skipped=glob.glob('*metadynLog_skipped')
    for f in glob.glob('sims*'):
        os.remove(f)
    for file in files_skipped:
        with open(file,'r') as ifile2:
            lines2=ifile2.readlines()
            split_file(file,lines2,args.parts)
    for j in range(args.parts):
        files_skipped_split=glob.glob(f'*metadynLog_skipped_{j}')
        files_split=glob.glob(f'*metadynLog_{j}')
        for file in files_split:
            #kappa,at=extract_params(file,args.pos,args.sep,args.restart)
            kappa,at=extract_params_restart(file,args.pos,args.sep,j,args.fac)
            #if at < float(args.start) or at > float(args.end):
            #    pass
                #print(f"Skipping {at} because outside boundaries")
            #else:
            with open(f'sims.txt','a') as ofile:
                    ofile.write(f'{file} {at} {kappa}\n')
        for file in files_skipped_split:
            #kappa,at=extract_params(file,args.pos,args.sep,args.restart)
            kappa,at=extract_params_restart(file,args.pos,args.sep,j,args.fac)
            #if at < float(args.start) or at > float(args.end):
            #    print(f"Skipping {at} because outside boundaries")
            #else:
            with open(f'sims_{j}.txt','a') as ofile:
                ofile.write(f'{file} {at} {kappa} 50 \n')
            with open(f'sims_equil.txt','a') as ofile:
                ofile.write(f'{file} {at} {kappa} 50 \n')
    do_wham(args.parts,f"{minwham:.3f}",f"{maxwham:.3f}",str(nbins),args.temp,args.tol,args.periodic)


if __name__=="__main__":
    main()
            
                
        
    #split_file(ifile,lines,5)
    


    
    
    
   
   
   