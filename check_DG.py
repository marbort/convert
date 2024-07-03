import sys
import argparse

parser = argparse.ArgumentParser(
                    prog='Check',
                    description='Calculate DG difference in kcalmol',
                    epilog='Energy is free')
parser.add_argument('--reacts', type=str, help='Reactants',nargs='+')
parser.add_argument('--freerads', type=str, help='Free rads',nargs='+')
parser.add_argument('--solvents', type=str, help='Solvents',nargs='+')
args = parser.parse_args()

data_react={}
data_freerad={}
for react in args.reacts:
    data[react]={'Neutral':{},'Radical':{}}
    for solvent in args.solvents:
            data[react][rad][solvent]

with open('SI.txt',r) as ifile:
    lines=ifile.readlines()

        