import sys
import glob
elm_dict={
    '1':'C',
    '2':'O',
    '3':'H',
    '4':'MG',
    '5':'CL'
}

def add_elm_pdb(file):
    with open(file,'r') as ifile:
        lines=ifile.readlines()
    newlines=[]
    for line in lines:
        if line[0:4] == "ATOM":
            elm=elm_dict[line[13]]
            letters=len(elm)
            newlines.append(line[:(-1*letters)-1]+elm)
        else:
            newlines.append(line)
            
    with open(file.replace('.pdb','_elm.pdb'),'w') as ofile:
        for line in newlines:
            ofile.write(f"{line}\n")


def main():
    files=glob.glob(sys.argv[1])
    for file in files:
        add_elm_pdb(file)


if __name__=="__main__":
    main()


            
        
    