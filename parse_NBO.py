import argparse
import os
import numpy
import pandas as pd
from operator import itemgetter

def extract_NBO_contributions(input):
    with open(input,'r') as ifile:
        lines=ifile.readlines()
    start=False
    start2=False
    start3=False
    NBOS_ok=False
    nbos=[]
    neda=[]
    for line in lines:
        if "Normal termination" in line:
            NBOS_ok=True
        if "NATURAL BOND ORBITALS (Summary):" in line:
            start=False
        if "Dipole moments:" in line:
            start2=False
        if "within unit  1" in line:
                if NBOS_ok:
                    pass
                else:
                    start=True
        if start:
            #print(line.split('.')[0])
            try: 
                int(line.split('.')[0])
                nbos.append([line[:30],line[30:59],float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])])
            except:
                pass
        if start2:
            if "Electrical" in line:
                start3=True
            if start3:
                try:
                    neda.append([line[0:23],float(line[25:37]),"","","",""])
                except:
                    neda.append([line[0:23],line[25:37],"","","",""])
            else:
                neda.append([line[0:15],line[15:34],line[34:53],line[53:60],line[60:71],line[72:79]])
        if "Natural Energy Decomposition Analysis (Summary)" in line:
            start2=True
    return(nbos,neda)
       

def convert_to_pandas(data,threshold,extract,name):
    header_nbos=["Donor", "Acceptor","E(2)", "E(NL)-E(L)", "F(L,NL)"]
    
    df_nbos=pd.DataFrame(data[0])
    df_neda=pd.DataFrame(data[1])
    df_nbos.columns=header_nbos
    writer = pd.ExcelWriter(f"{name}_NBOS.xlsx")   
    if extract:
        for k in range(0, len(extract),2):
            new=[]
            for row,i in enumerate(df_nbos["Donor"]):
                if f"{extract[k]}{int(extract[k+1]):3d}" in i:
                    new.append(df_nbos.loc[row])
            new_df=pd.DataFrame(new)

            try:
                new_df=new_df[new_df["E(2)"] > threshold]
                new_df.sort_values("E(2)",inplace=True,ascending=False)
                new_df.to_excel(writer,sheet_name=f"{extract[k]}{int(extract[k+1]):3d}")
            except:
                print(f"No contribution of atom {extract[0]}{extract[1]} found in NBO")
    df_nbos.sort_values("E(2)",inplace=True,ascending=False)
    df_filtered=df_nbos[df_nbos["E(2)"] > threshold]
    
    df_nbos.to_excel(writer,sheet_name="Full")
    df_filtered.to_excel(writer,sheet_name=f"Filtered {threshold}")
    df_neda.to_excel(writer,sheet_name="NEDA")
    writer.close()
    return(df_nbos,df_filtered,df_neda)


def main():
    parser = argparse.ArgumentParser(description='Plot data')
    parser.add_argument('--input', dest='input', 
                        type=str, help='input data')
    parser.add_argument('--output', dest='output', 
                        type=str, help='output data')
    parser.add_argument('--threshold', dest='threshold', 
                        type=float, help='E(2) threshold for printing')
    parser.add_argument('--extract', dest='extract', 
                        type=str, help='Extract donor atom',nargs='+')
    
    
    args = parser.parse_args()
    
    data=extract_NBO_contributions(args.input)

    df,df_filtered,df_neda=convert_to_pandas(data,args.threshold,args.extract,args.output)
   


if __name__=="__main__":
    main()
    
        
            