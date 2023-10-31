import json
import glob
import re
import sys
import logging
import os.path


def parse_tech(file,type):
    with open(file, 'r') as file:
        lines=file.read()
    if type == 'Gaussian':    
        match=re.findall("\\\\#p.*?\\\\",lines,re.DOTALL)
    elif type == 'ORCA':
        match=[re.sub("\| *\d*>","",x).encode('unicode_escape') for x in re.findall("(?<=1>).*?(?=\*xyz)",lines,re.DOTALL)]
    return(match,len(match))
    



                
            
    
    
    
    
