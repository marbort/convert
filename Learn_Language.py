import os
import argparse
import numpy as np



def give_words(language,number):
    try:
        with open(f'{language}_remaining_words.txt','r') as ifile:
            words=[x[0] for x in ifile.readlines().split()]
    except:
        with open(f'{language}_all_words.txt','r') as ifile:
            words=[x[0] for x in ifile.readlines().split()]
    selected=np.random.random_integers(0,len(words),number)
    
    
    
    
    with open(f'{language}_learned_words.txt','a') as ofile:
            for i in selected:
                ofile.write(words[i])
    with open(f'{language}_remaining_words.txt','a') as ofile2:
        for j,word in enumerate(words):
            if j in selected:
                pass
            else:
                ofile2.write(word)
            




give_words(10)
            
    
    
        