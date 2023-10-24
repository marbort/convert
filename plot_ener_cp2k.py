import matplotlib.pyplot as plt
import argparse
import re


parser = argparse.ArgumentParser(description='Plot MD ener data CP2K')


parser.add_argument('--input' , dest='input',help='ener file')
parser.add_argument('--skip' , dest='skip',default=0,type=int,help='skip n first y columns')
args = parser.parse_args()

with open(args.input,'r') as ifile:
    lines=ifile.readlines()
data={}
header=re.split(r'\s{2,}', lines[0])
for i in range(1,len(header)):
    data[header[i]]=[]
#print(data)
for line in lines[1:]:
    for i,item in enumerate(data):
        data[item].append(float(line.split()[i]))

fig=plt.figure(figsize=(16,9),dpi=150)
plt.suptitle(args.input)
for i in range(4):
    plt.subplot(2,2,i+1)
    plt.plot(data[list(data.keys())[0]],data[list(data.keys())[i+1+args.skip]])
    plt.xlabel(list(data.keys())[0])
    plt.ylabel(list(data.keys())[i+1+args.skip])
    #plt.title(list(data.keys())[i+2])
plt.savefig('data_plot.png',format='png')
    
