with open('colvar','r') as ifile:
    lines=ifile.readlines()
newlines=[]
newlines.append(lines[0])
print(len(lines[0].split()))    
print(len(lines[1].split()))    
for line in lines[1:]:
    if len(line.split()) == len(lines[0].split())-2:
        newlines.append(f'{line.strip()}\n')
    else:
        newlines.append(f'{line.strip()} 0.0\n')
with open('colvar_ok','w') as ofile:
    ofile.writelines(newlines)
    
