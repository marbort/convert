import os


def convert_tbut_et(file):
    idx=[7,8,9,15,16,17,19,20,21,49,50,51]
    idx_h=[6,14,18,48]
    converted=[]
    with open(file,'r') as ifile:
        lines=ifile.readlines()
    for i,line in enumerate(lines):
        if i in idx:
            continue
        elif i in idx_h:
            converted.append(line.replace("C","H"))
        else:
            converted.append(line)
    return(converted)
            
        


for file in os.listdir('./'):
    if 'react' in file:
        #print(file)
        converted=convert_tbut_et(file)
        #print(file.replace('tBu','Et'))
        with open(file.replace('tBu','Et'),'w') as ofile:
            #print(ofile)
            for line in converted:
               ofile.write(line)