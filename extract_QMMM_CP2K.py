import argparse
import itertools
import json
 
def intervals_extract(iterable):
     
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
    lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]


parser = argparse.ArgumentParser(description='Plot data')


parser.add_argument('--input', dest='input', default=None,
                    type=str, help='PDB file containing the whole system structure')
parser.add_argument('--qm', dest='qm', type=str,
                    help='PDB file containing the QM cluster')
args = parser.parse_args()
cluster=[]
natoms=0

with open(args.qm,'r') as clusterfile:
    lines=clusterfile.readlines()
for line in lines[1:-1]:
    #cluster.append("".join(["{:8.3f}".format(float(line.split()[x])) for x in range(6,9)]))  
    if line[22:26] not in cluster:
        cluster.append(f"{line[22:26]:4s}")
    natoms=int(line.split()[1])

with open(args.input,'r') as ifile:
    lines=ifile.readlines()
print(cluster)
res=[1,
5,
294,
1923,
2311,
2572,
3236,
3298,
3435,
3512,
3622,
3792,
4170,
4405,
5797]

##################
#######NEED TWEAKING########
######################


QMmm_atoms={}
elms=["H","C","O","N","MG","CL"]
for elm in elms:
    QMmm_atoms[elm]=[]
for i,line in enumerate(lines):
        if f"{line[22:26]:4s}" in cluster:
            QMmm_atoms[line.split()[-1].rstrip()].append((line.split()[1],i))
            #break
with open ('QMM_atoms.json','w') as jfile:
    json.dump(QMmm_atoms,jfile,indent=2)

atom_idx=[]
for elm in QMmm_atoms:
    for idx in QMmm_atoms[elm]:
        if int(idx[0]) == int(idx[1]):
            atom_idx.append(int(idx[0])-1)
        else:
            print(f"Error in atom {idx[0]} found at line {idx[1]}")
atom_idx.sort()
intervals=list(intervals_extract(atom_idx))
with open('vmd_index.dat','w') as ofile:
    for interval in intervals:
        try:
            ofile.write(f"{interval[0]} to {interval[1]} ")
        except:
            ofile.write(f"{interval[0]} ")

extracted=sum([len(QMmm_atoms[x]) for x in QMmm_atoms])
if int(extracted) != int(natoms):
    print('Error in atom numbers: input {} extracted {}'.format(natoms,extracted))
    print(QMmm_atoms)

"""


index={}
for elm in elms:
    index[elm]=[]
#print(index)
for i in res:
    for j in lines:
        try:
            #print(j[22:26].lstrip())
            if str(i) == j[22:26].lstrip():
                QMMM_atoms.append(j)
        except:
            continue


for line in QMMM_atoms:
    index[line.split()[-1]].append(line.split()[1])
tot_len=[len(index[x]) for x in index]  
print(index)
print(sum(tot_len))
"""
with open('QMMM_atoms.dat','w') as ofile:
    for i in QMmm_atoms:
        ofile.write("&QM_KIND {}\n ".format(i))
        ofile.write("MM_INDEX ")
        for j in QMmm_atoms[i]:
            ofile.write(f"{j[0]} ")
        ofile.write("\n")
        ofile.write("&END QM_KIND\n")

           
        

    