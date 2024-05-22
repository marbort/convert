import numpy as np
import sys

#type_map={"C":1,"O":2,"H":3,"Mg":6,"Cl":4,"Li":5}
type_map={"C":1,"O":2,"H":3,"Mg":4,"Cl":5,"Li":6}


elm,x,y,z=np.loadtxt(sys.argv[1],unpack=True,skiprows=2,dtype=str)
types=[type_map[x] for x in elm]

with open('conf_data_temp','w') as ofile:
    for i,line in enumerate(types):
        ofile.write("{:8d}\t{:8d}\t{:10s}\t{:10s}\t{:10s}\n".format(i+1,line,x[i],y[i],z[i]))





print('ciao')
