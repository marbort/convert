import argparse
import os
import numpy
import pandas as pd
from operator import itemgetter

def extract_EDA(input):
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    start=False
    start2=False
    EDA=[]
    NOCV=[]
    idx=0
    for line in lines:
        if "Correction terms (incorporated" in line:
            start=False
        if "SFO decomposition of Delta rho k (major contributions)" in line:
            start2=False
        if "on the meaning of the various terms" in line:
            idx += 1
            if idx == 3:
                start=True
        if start:
            #print(line.split('.')[0])
            try: 
                float(line.split()[-1])
                EDA.append([line[:38],float(line.split()[-4]),float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])])
            except:
                pass
        if "Orbital Interaction Energy Contributions from each NOCV pair (in kcal/mol)" in line:
            start2=True
        if start2:
            try:
                float(line.split()[-1])
                NOCV.append([line.split()[0],float(line.split()[1])])
            except:
                pass
    NOCV.append(["Total Sum",sum([x[1] for x in NOCV])])
    return(EDA,NOCV)


def convert_to_pandas(data,writer):
    header_EDA=["Term", "hartree", "eV", "kcal/mol", "kJ/mol"]
    header_NOCV=["Pair","Interaction"]
    df_EDA=pd.DataFrame(data[0])
    df_EDA.columns=header_EDA
    df_NOCV=pd.DataFrame(data[1])
    df_NOCV.columns=header_NOCV
    df_EDA.to_excel(writer,sheet_name="EDA")
    df_NOCV.to_excel(writer,sheet_name="NOCV")
    


def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--input', dest='input', 
                        type=str, help='input data')
    parser.add_argument('--output', dest='output', 
                        type=str, help='input data')
    parser.add_argument('--threshold', dest='threshold', 
                        type=float, help='E(2) threshold for printing')
    parser.add_argument('--extract', dest='extract', 
                        type=str, help='Extract donor atom',nargs='+')
    
    args = parser.parse_args()
    
    data=extract_EDA(args.input)
    writer = pd.ExcelWriter(f'{args.output}_ESTNOCV.xlsx') 
    convert_to_pandas(data,writer)
    writer.close()
   


if __name__=="__main__":
    main()
    
        
            