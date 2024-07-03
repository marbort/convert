import numpy as np
import sys


def extract_data(file):
    data={}
    with open(file,'r') as ifile:
        lines=ifile.readlines()
        data['file']=ifile.name
    for line in lines:
        if line[0] == '#':
            continue
        if 'title' in line:
            data['title']=" ".join(line.split()[2:len(line.split())])
        if 'xaxis  label' in line:
            data['xlabel']=" ".join(line.split()[3:len(line.split())])
        if 'yaxis  label' in line:
            data['ylabel']=" ".join(line.split()[3:len(line.split())])
        if 's0 legend' in line:
            data['label']=" ".join(line.split()[3:len(line.split())])
    print("Found data with legend: {}".format(data['label']))
    data['x']=[float(x.split()[0]) for x in lines if "#" not in x if "@" not in x]
    data['y']=[float(x.split()[1]) for x in lines if "#" not in x if "@" not in x]
    return(data)





data=extract_data(sys.argv[1])

avg=np.average(data['y'])
std=np.std(data['y'])


print(avg,std)