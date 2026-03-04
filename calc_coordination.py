import numpy as np
import argparse

def calc_coordination(dist, cutoff=3.0,exponent=6):
    """
    Calculate the coordination number based on distance and cutoff.
    """
    a=1-(dist/cutoff)**exponent
    b=(1-(1/cutoff)**(exponent*2))
    print(a,b)
    return (1-(dist/cutoff)**exponent)/(1-(dist/cutoff)**(exponent*2))


def main():
    parser = argparse.ArgumentParser(description='Plot a slice of free energy surface.')
    parser.add_argument('bonds', type=str, nargs='+',
                        help='List of bond files to process')
    parser.add_argument('--cutoff', type=float,  default=3.0)
    
    parser.add_argument('--exponent', type=float, default=6.0,
                        help='Exponent for coordination number calculation')
    
    args= parser.parse_args()
    
    CN_arr=[calc_coordination(float(i), cutoff=parser.parse_args().cutoff, exponent=parser.parse_args().exponent) for i in args.bonds]
    print(CN_arr)
    print(f"Coordination numbers: {sum(CN_arr)}")
    
if __name__ == "__main__":
    main()