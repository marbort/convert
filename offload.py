import subprocess
import os
import timefoldercp
import json

def offload_calculation(cluster,prog):
    slurm_ids=[]
    input_basename=os.path.basename(input_dict['offload'][prog][cluster]['input'])
    data_basename=os.path.basename(input_dict['offload'][prog]['data'])
    #copy data
    subprocess.call(["scp","-i",input_dict['offload'][cluster]['ssh-key'],"-r",input_dict['offload'][prog]['data'],
                     input_dict['offload'][cluster]['user']+"@"+input_dict['offload'][cluster]['ip']+":"+input_dict['offload'][cluster]['folder']])
    #copy SLURM input
    subprocess.call(["scp","-i",input_dict['offload'][cluster]['ssh-key'],"-r",input_dict['offload'][prog][cluster]['input'],
                     input_dict['offload'][cluster]['user']+"@"+input_dict['offload'][cluster]['ip']+":"+os.path.join(input_dict['offload'][cluster]['folder'],data_basename)])
    #run SLURM
    job_id = subprocess.check_output(["ssh","-i",input_dict['offload'][cluster]['ssh-key'],
                     input_dict['offload'][cluster]['user']+"@"+input_dict['offload'][cluster]['ip'],
                     "cd",os.path.join(input_dict['offload'][cluster]['folder'],data_basename),";", "sbatch",input_basename]).decode('utf-8').split()[-1]
    slurm_ids.append(job_id)

    #check SLURM
    while True:
            time.sleep(60)
            all_finished = True
            for job_id in slurm_ids:
                # ignore the first line of squeue output
                try:
                    if subprocess.check_output(["ssh","-i",input_dict['offload'][cluster]['ssh-key'],
                                                input_dict['offload'][cluster]['user']+"@"+input_dict['offload'][cluster]['ip'],"squeue", "-j", job_id]).decode('utf-8').splitlines()[1:]:
                        all_finished = False
                        break
                except:
                    pass
            if all_finished:
                subprocess.call(["scp","-i",input_dict['offload'][cluster]['ssh-key'], "-r",
                     input_dict['offload'][cluster]['user']+"@"+input_dict['offload'][cluster]['ip']+":"+os.path.join(input_dict['offload'][cluster]['folder'],data_basename),
                                                                                                             "."])
                subprocess.call(["ssh","-i",input_dict['offload'][cluster]['ssh-key'], 
                     input_dict['offload'][cluster]['user']+"@"+input_dict['offload'][cluster]['ip'],
                     "rm","-rf",os.path.join(input_dict['offload'][cluster]['folder'],data_basename)])                                                                                                             "."])
                break
    
    

with open("input_dict.json",'r') as ifile:
    input_dict=json.load(ifile)

offload_calculation('lumi','cp2k')
    

