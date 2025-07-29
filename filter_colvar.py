import sys
import numpy as np
import shutil
import argparse
import pandas as pd


parser = argparse.ArgumentParser(description='Filter colvar file')
parser.add_argument('--colvar', type=str, help='Colvar file to filter', default='colvar')
parser.add_argument('--col', type=str, help='Columns to filter', default="CV1")
parser.add_argument('--valrange', type=float, help='Value to filter', nargs=2, default=[0,1])
parser.add_argument('--stride', type=float, help='Stride to filter', default=1)
parser.add_argument('--keep', help='Keep additional columns', default=[] ,nargs='+')

args=parser.parse_args()
filename = args.colvar
with open(filename, 'r') as f:
    fields = f.readline()
    colnames = fields.split()[2:]
print(colnames)
if args.keep[0] == "all":
    data=pd.read_table(filename,dtype=float,sep='\s+',comment='#',header=None,names=colnames)
    data=data.iloc[::int(args.stride),:]
    colnum=colnames.index(args.col)
    data_filtered=data[(data.iloc[:,colnum] >= args.valrange[0]) & (data.iloc[:,colnum] <= args.valrange[1])]
    data_filtered.to_csv(f"{filename}_{args.col}_{args.valrange[0]}to{args.valrange[1]}", sep=' ', header=False, index=False)
    with open(f"{filename}_{args.col}_{args.valrange[0]}to{args.valrange[1]}", 'r+') as file:
        content = file.read()
        file.seek(0, 0)
        file.write(fields.rstrip('\r\n') + '\n' + content)
else:
    data=pd.read_table(filename,dtype=float,sep='\s+',comment='#',header=None,names=colnames,usecols=["time",args.col]+args.keep)
    data=data.iloc[::int(args.stride),:]
    data.reset_index(drop=True, inplace=True)
    data_filtered=data[(data[args.col] >= args.valrange[0]) & (data[args.col] <= args.valrange[1])]
    data_filtered.to_csv(f"{filename}_{args.col}_{args.valrange[0]}to{args.valrange[1]}", sep=' ', header=True, index=True)

"""
shutil.copy(sys.argv[1],sys.argv[1]+'.bak')

with open(sys.argv[1], 'r') as f:
    lines = f.readlines()

col=lines[0].split().index("CV2")-2
print(col)
newlines    = []
for line in lines:
    if line.startswith('#'):
        newlines.append(line)
    elif float(line.split()[col]) < 1.5:
        pass
    else:
        newlines.append(line)
with open(sys.argv[1], 'w') as f:
    for line in newlines:
        f.write(line)
"""
    
