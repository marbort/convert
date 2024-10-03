import sys
import glob

paths=glob.glob(sys.argv[1])
#print(paths)
for file in paths:

    with open(file,'r') as ifile:
        lines=ifile.readlines()

    for i,line in enumerate(lines):
        if "CHARGE" in line:
            chrg_input=line.split()[-1]
        if "&QM_KIND CL" in line:
            CL_numb=len(lines[i+1].split()[1:])
        if "&QM_KIND N" in line:
            N_numb=len(lines[i+1].split()[1:])

    tot_charge_qm_atoms=2-CL_numb+N_numb
    try:
        assert int(chrg_input) == tot_charge_qm_atoms  
    except:
        print(f"Charge in input {file} should be {tot_charge_qm_atoms}")
        
