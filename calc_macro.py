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

def calculate_volume(density,MW,molecules):
    Navog=6.022e23
    nmols=molecules/Navog
    VmL=MW*nmols/density
    Vang=VmL*1e24
    return(Vang,(Vang)**(1/3),round(Vang**(1/3)-1,0))


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

if sys.argv[1] == "volume":
    vml,side,packside=calculate_volume(float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]))
    print("=========================================")
    print(f"Density (g/mL)= {sys.argv[2]} MW= {sys.argv[3]} Molecules= {sys.argv[4]}")
    print(f"Volume (Å^3)= {vml:.4f}")
    print(f"Side (Å)= {side:.2f}")
    print(f"packmol Side (Å)= {packside:.2f}")
    print("=========================================")