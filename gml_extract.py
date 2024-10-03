import glob
import pandas as pd

files=glob.glob('*.gml')


dict={"ID":[],"atoms":[],"metal":[]}


for file in files:
    with open(file) as fin:
        fin.seek(0)
        data = fin.read(250)
        dict["ID"].append(data[45:51])
        dict["atoms"].append(data[65:67]) #atoms
        dict["metal"].append(data[114:116]) #metal
print(dict)
df=pd.DataFrame.from_dict(dict)

print(df)
        