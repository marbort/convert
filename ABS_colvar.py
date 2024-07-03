import pandas as pd
import sys

with open(sys.argv[1],'r') as ifile:
    data=pd.read_csv(ifile,sep=" ",header=0)
data['CV1']=abs(data['CV1'])
#data=data.loc[data['cvMg1THF'] <= 3.5]
pd.DataFrame.to_csv(data,"colvar_abs",sep=" ",index=False)


print(data) 

