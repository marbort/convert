from inspect import CO_ASYNC_GENERATOR
import numpy as np
import sys
import argparse

#parser = argparse.ArgumentParser(description='Exctract xyz coordinates')
#parser.add_argument()
#parser.add_argument('--type', dest='tipo', default='macro',  help='Calculation type: macro or nmols')
#args = parser.parse_args()

def calculate_macroscopic(nmols,boxx,boxy,boxz,MW):
    Navog=6.022e23
    amu=1.66054e-24
    Vang=boxx*boxy*boxz
    VmL=Vang*1e-24
    VL= VmL/1000
    Conc=nmols/Navog/VL
    mass=MW*nmols
    Density=mass*amu/VmL

    return(Conc, Density,VmL,mass)

def calculate_nmols(boxx,boxy,boxz,d,MW):
    Navog=6.022e23
    amu=1.66054e-24
    Vang=boxx*boxy*boxz
    VmL=Vang*1e-24
    nmols=VmL*d*Navog/MW
    Density=MW*nmols*amu/VmL

    return(nmols,VmL)

def calculate_macroscopic_mass(mass,boxx,boxy,boxz):
    amu=1.66054e-24
    Vang=boxx*boxy*boxz
    VmL=Vang*1e-24
    VL= VmL/1000
    Density=mass*amu/VmL
    return(Density,VmL)


if sys.argv[1] == "macro": 
    Conc,Density,vml,mass=calculate_macroscopic(float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]) \
                    ,float(sys.argv[5]),float(sys.argv[6]))

    print("""
    Concentration= {:.3f} M
    Density= {:.3f} g/mL
    """.format(Conc,Density))

if sys.argv[1] == "nmols":
    nmols=calculate_nmols(float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),\
        float(sys.argv[5]),float(sys.argv[6]))
    print(f"nmols= {nmols[0]:.0f} ")
    print(f"V (mL)= {nmols[1]:.2e} ")

if sys.argv[1] == "macro_mass":
    Density,vml=calculate_macroscopic_mass(float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]),\
        float(sys.argv[5]))
    print("density= {},{} ".format(Density,vml))