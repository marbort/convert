import glob

def load_data(inputs):
    lines=[]
    for input in inputs:
        print(f"Extracting {input}")
        with open(input,'r') as ifile:
            tmp=ifile.readlines()
            for line in tmp:
                lines.append(line)
    return(lines)

def write_data(lines,offset,outfile):
    natoms=int(lines[3])
    header=9
    totframes=len(lines)//(natoms+header)
    newframes=totframes//offset
    newlines=[]
    for i in range(newframes):
        for j in range(natoms+header):
            newlines.append(lines[i*(natoms+header)*offset+j])
    with open(outfile,'w') as ofile:
        for line in newlines:
            ofile.write(line)


inputs=glob.glob("window_*/lmp.lammpstrj")
lines=load_data(inputs)
write_data(lines,10,"newlmp.lammpstrj")
            
    








