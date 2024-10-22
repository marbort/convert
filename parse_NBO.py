import argparse
import os
import numpy
import pandas as pd
from operator import itemgetter

def extract_NBO_contributions(input):
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    start=False
    nbos=[]
    for line in lines:
        if "NATURAL BOND ORBITALS (Summary):" in line:
            start=False
        if "within unit  1" in line:
            start=True
        if start:
            #print(line.split('.')[0])
            try: 
                int(line.split('.')[0])
                nbos.append([line[:30],line[30:59],float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])])
            except:
                pass
    return(nbos)

def convert_to_pandas(nbos,threshold,extract):
    header=["Donor", "Acceptor","E(2)", "E(NL)-E(L)", "F(L,NL)"]
    df=pd.DataFrame(nbos)
    df.columns=header
    writer = pd.ExcelWriter('NBOS.xlsx')   
    if extract:
        for k in range(0, len(extract),2):
            new_df=pd.DataFrame()
            for row,i in enumerate(df["Donor"]):
                if f"{extract[k]}{int(extract[k+1]):3d}" in i:
                    new_df=new_df.append(df.loc[row])
            new_df=new_df[new_df["E(2)"] > threshold]
            new_df.sort_values("E(2)",inplace=True,ascending=False)
            new_df.to_excel(writer,sheet_name=f"{extract[k]}{int(extract[k+1]):3d}")
    df.sort_values("E(2)",inplace=True,ascending=False)
    df_filtered=df[df["E(2)"] > threshold]
    
    df.to_excel(writer,sheet_name="Full")
    df_filtered.to_excel(writer,sheet_name=f"Filtered {threshold}")
    writer.close()
    return(df,df_filtered)


def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--input', dest='input', 
                        type=str, help='input data')
    parser.add_argument('--threshold', dest='threshold', 
                        type=float, help='E(2) threshold for printing')
    parser.add_argument('--extract', dest='extract', 
                        type=str, help='Extract donor atom',nargs='+')
    
    args = parser.parse_args()
    
    nbos=extract_NBO_contributions(args.input)
    df,df_filtered=convert_to_pandas(nbos,args.threshold,args.extract)
   


if __name__=="__main__":
    main()
    
        
            