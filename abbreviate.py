import os
import json
import regex as re
import sys

def create_database(list,outpath):
    database={}
    with open(list,'r') as ifile:
        lines=ifile.readlines()
        for line in lines:
            abbrv=[x+"." for x in line.split('=')[1].rstrip().split()]
            abbrvdot=' '.join(abbrv)
            dat_name=re.sub('[^a-zA-Z0-9 \n\.]', '', line.split('=')[0].rstrip().lower())
            database[dat_name]=abbrvdot
    with open(os.path.join(outpath,'abbreviations.json'),'w') as ofile:
        json.dump(database,ofile,indent=1)

def abbreviate_bib(bibfile,database):
    with open(database,'r') as datafile:
        abbrv=json.load(datafile)
    with open(bibfile,'r') as ifile:
        lines=ifile.readlines()
        for i,line in enumerate(lines):
            if "journal" in line:
                jname=re.findall('\{(.*?)\}',line)
                try:   
                    tmp=line.replace(jname[0],abbrv[jname[0].lower().rstrip()])
                    lines[i]=tmp
                except:
                    try:
                        tmp=line.replace(jname[0],abbrv[jname[0][4:].lower().rstrip()])
                        lines[i]=tmp
                    except:
                        print("Journal {} not found".format(jname))
        with open('bib_abbrv.bib','w') as ofile:
            for line in lines:
                ofile.write(line)


#create_database('/home/marco/Downloads/jabref_wos_abbrev.txt','/home/marco')
abbreviate_bib(sys.argv[1],'/home/marco/abbreviations.json')
                    
        
    
            