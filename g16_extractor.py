import cclib
import sys
import os
import json

"""
EXAMPLE JSON FILE 

{
  "files":[
             "mp2.log",
             "b2plyp.log"
             ],

  "properties":[
             "natom",
             "scfenergies",
             "mpenergies"
            ]

 }

"""


input = sys.argv[1]
lines=[]
with open(input,'r') as ifile:
    parms=json.load(ifile)
for file in parms['files']:
    parser=0
    data=0
    parser = cclib.io.ccopen(file)
    data = parser.parse()
    tmp=[]
    for prop in parms['properties']:
        try:
            tmp.append(getattr(data,prop)[-1])
        except:
            if getattr(data,prop):
                tmp.append(getattr(data,prop))
            else:
                tmp.append(0)
    lines.append([file]+tmp)
with open ('results.dat', 'w') as ofile:
    ofile.write(",".join(["file"]+parms['properties']))
    ofile.write("\n")
    for line in lines:
        ofile.write("{}\n".format(",".join([str(x) for x in line])))
    
print(lines)
    
    
