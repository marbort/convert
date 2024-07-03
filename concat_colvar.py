import argparse
import numpy as np
import pandas as pd
import sys
import plumed as pld
 


parser = argparse.ArgumentParser(description='Format colvar file for use with WHAM')

parser.add_argument('--input',dest='input',type=str,default='colvar',help='files to concatenate',nargs='+')
parser.add_argument('--output',dest='output',default='colvar_full',type=str,help='concatenated file')


args = parser.parse_args()

files=[]
for input in args.input:
    files.append(pld.read_as_pandas(input,enable_constants=False))
concat=pd.concat(files,join='outer')
print(concat.shape)

concat_filled = concat.fillna(0)

with open(args.output,'w') as ofile:
    ofile.write(f"#! FIELDS {' '.join(concat_filled.columns.values)}\n")<
    concat_filled.to_csv(ofile,sep=" ",float_format='%6.4f',index=False,header=False)
"""
with open(args.output,'w') as ofile:
    ofile.write("#! FIELDS")
    for n in concat_filled.columns:
            ofile.write(" "+str(n))
    ofile.write("\n")
    #for i in range(2):
    for i in range(concat_filled.shape[0]):
        for j in concat_filled.columns:
     #       print(concat_filled[j][i].values[0])
            #try:
            #    ofile.write(f" {concat_filled[j][i].values[0]:10.4f}")
            #except:
                ofile.write(f" {concat_filled[j][i]:10.4f}")
        ofile.write("\n")

"""

#print(concat_filled)
#print(files[0])
#pld.write_pandas(concat_filled,f"{args.output}")
#pld.write_pandas(files[1],f"{args.output}")
    