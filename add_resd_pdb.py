import  numpy
import sys



#division=[(0,12,12,"IMC"),(12,2560,13,"THF")] Me
#division=[(0,13,13,"IMC"),(13,40,13,"THF"),(39,44,5,"IMC"),(44,2566,13,"THF")] #Et
division=[(0,24,24,"IMC"),(24,2572,13,"THF")] #iPr

with open(sys.argv[1],'r') as ifile:
    lines=ifile.readlines()

newlines=[]
idx=1
for div in range(len(division)):  
    for j in range((division[div][1]-division[div][0])//division[div][2]):
        for k in range(j*division[div][2],j*division[div][2]+division[div][2]):
            tmpline=lines[division[div][0]+k+2].replace("UNK",division[div][3])
            tmpline=tmpline[:22]+f"{idx:4d}"+tmpline[27:]
            newlines.append(tmpline)
        idx+=1

with open('test.pdb','w') as ofile:
    for line in newlines:
        ofile.write(line)


            
            
            