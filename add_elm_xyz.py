import argparse
import dpdata as dp

parser = argparse.ArgumentParser(description="Plot data")
parser.add_argument("--input", dest="input", type=str, help="input data")
parser.add_argument(
    "--map", dest="map", type=str, help="Mapping from types to elements", default="TurboMarinella"
    )

args = parser.parse_args()

maps= {"TurboMarinella":{1:"C",2:"O",3:"H",4:"Cl",5:"Li",6:"Mg"}}


data=dp.System(args.input,'xyz')

for i,name in enumerate(data['atom_names']):
    data['atom_names'][i]=maps[args.map][int(name)]

data.to_xyz('new.xyz')
    
        
        


