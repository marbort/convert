import os
import argparse
import numpy as np
from deep_translator import GoogleTranslator



def give_words(language,number):
    try:
        with open(f'{language}_remaining_words.txt','r') as ifile:
            words=[x[0] for x in ifile.readlines().split()]
    except:
        with open(f'{language}_all_words.txt','r') as ifile:
            words=[x.split()[0] for x in ifile.readlines()]
    selected=np.random.randint(0,len(words),number)
    translated=[GoogleTranslator(source=f'{language.lower()}', target='en').translate(text=words[x]).lower() for x in selected]
    translated_all=[GoogleTranslator(source=f'{language.lower()}', target='en').translate(text=x).lower() for x in words]
    print(selected)
    print(translated)
    
    
    with open(f'{language}_all_words.txt','w',encoding="utf-8") as allfile:
        for i,word in enumerate(words):
            allfile.write(f"{word}\t{translated_all[i]}")
    
    with open(f'{language}_current_words.txt','w',encoding="utf-8") as ofile:
            for i,idx in enumerate(selected):
                ofile.write(f'{words[idx]}\t{translated[i]}\n')
    with open(f'{language}_current_words_exercise.txt','w',encoding="utf-8") as ofile:
            for i,idx in enumerate(selected):
                ofile.write(f'{words[idx]}\n')
    
    with open(f'{language}_learned_words.txt','a') as ofile2:
            for i,idx in enumerate(selected):
                ofile2.write(f'{words[idx]}\t{translated[i]}')
    with open(f'{language}_remaining_words.txt','a') as ofile3:
        for j,word in enumerate(words):
            if j in selected:
                pass
            else:
                ofile3.write(word)

def check_words(language,exercise):
    with open(exercise,'r') as ifile:
        lines=ifile.readlines()
        original=[x.split()[0].lower() for x in lines]
        translated=[" ".join(x.split()[1:]).lower() for x in lines]
    with open(f'{language}_current_words.txt','r') as checkfile:
        check=[" ".join(x.split()[1:]).lower() for x in checkfile.readlines()]
    print(translated,original)
    for i in range(len(translated)):
        if translated[i]==check[i]:
            pass
        else:
            print(f"Translation of {original[i].upper()} is not {translated[i].upper()} but {check[i].upper()}")

                  
        


def main():
    
    parser = argparse.ArgumentParser(description='Plot data')

    parser.add_argument('--language' , dest='language',default='italian', type=str,help='Language you want to learn')
    parser.add_argument('--words' , dest='words',default=10,type=int,help='Number of words to learn')
    
    
    args = parser.parse_args()




    give_words(args.language,args.words)
    if args.check:
        check_words(args.language,"Norwegian_current_words_exercise.txt")
                
    
if __name__=="__main__":
    main()
        