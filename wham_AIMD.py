
import argparse
import glob
import os
import subprocess
import sys

def skip_lines(file,lines,skip):
    with open(f"{file.name}_skipped",'w') as ofile:
        for line in lines[skip:]:
            ofile.write(line)

def extract_params(file,pos,sep,restart):
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
    return(float(kappa)*2,round(float(at)/0.5292,2))

def extract_params_restart(file,pos,sep,idx):
    at=file.split('-')[0].split(f'{sep}')[pos]
    cp2k_restart=file.replace(f"-COLVAR.metadynLog_skipped_{idx}",'-1.restart')
    kappa=0
    try:
        with open(cp2k_restart,'r') as cpinp:
            lines=cpinp.readlines()
        for line in lines:
            if 'K  ' in line:
                kappa=line.split()[-1]
        return(round(float(kappa)*2*2625.5002,1),round(float(at)/0.5292,2))
    except Exception as X:
        print(X)
        print(f"Skipping {at} restart. File not found")
        return(round(float(kappa)*2*2625.5002,1),round(float(at)/0.5292,2))
    
    

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
                   
def do_wham(parts,start,end,nbins,temp):
    subprocess.run(["/home/marco/WHAM/wham/wham/wham", start, end, nbins, "0.00001", temp, "0", "sims.txt", "out.dat_all"])
    subprocess.run(["/home/marco/WHAM/wham/wham/wham", start, end, nbins, "0.00001", temp, "0", "sims_equil.txt", "out.dat_all_equil", "50", "3153"])
    for i in range(parts):
       subprocess.run(["/home/marco/WHAM/wham/wham/wham", start, end, nbins, "0.00001", temp, "0", f"sims_{i}.txt", f"out.dat_all_{i}"])
    
    
    


def main():
    
    parser = argparse.ArgumentParser(description='Plot data')

    parser.add_argument('--parts' , dest='parts',default=5, type=int,help='Number of splits')
    parser.add_argument('--pos' , dest='pos',default=4,type=int,help='Restrain position in filename')
    parser.add_argument('--sep' , dest='sep',default='_',help='Separator to get restrain position')
    parser.add_argument('--restart' , dest='restart',action='store_true', default=False, help='Use CP2K restart for restraint parameters',)
    parser.add_argument('--skip' , dest='skip',default=1000, type=int,help='Number of lines to skip from colvar file')
    parser.add_argument('--start' , dest='start',type=str,default=0,help='Start of CV')
    parser.add_argument('--end' , dest='end',default="1",type=str,help='End of CV')
    parser.add_argument('--nbins' , dest='nbins', type=str, default="50",help='Number of WHAM bins')
    parser.add_argument('--temp' , dest='temp',type=str, default="300",help='Temperature for WHAM')
    parser.add_argument('--fac' , dest='fac',type=float, default=1.89,help='Conversion factor between at in filename and in colvar. Default: 1.89')
    
    
    args = parser.parse_args()
    
    files=glob.glob('*metadynLog')
    
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
            kappa,at=extract_params_restart(file,args.pos,args.sep,j)
            if at < float(args.start) or at > float(args.end):
                pass
                #print(f"Skipping {at} because outside boundaries")
            else:
                with open(f'sims.txt','a') as ofile:
                    ofile.write(f'{file} {at} {kappa}\n')
        for file in files_skipped_split:
            #kappa,at=extract_params(file,args.pos,args.sep,args.restart)
            kappa,at=extract_params_restart(file,args.pos,args.sep,j)
            if at < float(args.start) or at > float(args.end):
                print(f"Skipping {at} because outside boundaries")
            else:
                with open(f'sims_{j}.txt','a') as ofile:
                    ofile.write(f'{file} {at} {kappa} 200 \n')
                with open(f'sims_equil.txt','a') as ofile:
                    ofile.write(f'{file} {at} {kappa} 200\n')
    do_wham(args.parts,args.start,args.end,args.nbins,args.temp)


if __name__=="__main__":
    main()
            
                
        
    #split_file(ifile,lines,5)
    


    
    
    
   
   