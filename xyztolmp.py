import numpy as np
import sys

type_map={"Turbo":{"C":1,"O":2,"H":3,"Mg":4,"Cl":5,"Li":6},
          "Normal":{"C":1,"O":2,"H":3,"Mg":4,"Cl":5},
          "Hauser":{"C":1,"O":2,"H":3,"Cl":4,"Li":5,"Mg":6, "N":7},
          "Turbo_Marinella":{"C":1,"O":2,"H":3,"Cl":4,"Li":5,"Mg":6,"Br":7, "I":8}
}

#type_map={"C":1,"Cl":2,"H":3,"N":4,"O":5,"Mg":6}



elm,x,y,z=np.loadtxt(sys.argv[1],unpack=True,skiprows=2,dtype=str)
types=[type_map[sys.argv[2]][x] for x in elm]

with open('conf_data_temp','w') as ofile:
    for i,line in enumerate(types):
        ofile.write("{:8d} {:8d} {:10.6f} {:10.6f} {:10.6f}\n".format(i+1,line,float(x[i]),float(y[i]),float(z[i])))





print('ciao')
