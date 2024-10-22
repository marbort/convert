import argparse



def get_sigma(input,colvars):
    idx=[]
    sigmas=["sigma_"+x for x in colvars]
    sigma=[]
    with open(input,'r') as ifile:
        lines=ifile.readlines()
        header=lines[0]
        last=lines[-3:-1]
    for i,col in enumerate(header.split()):
        if col in sigmas:
            idx.append(i-2)
    for ind in idx:
        try:
            sigma.append(last[-1].split()[ind])
        except:
            print("KERNEL file missing data")
    return(sigma)
            
        
        
    
    
def main():
    
    parser = argparse.ArgumentParser(description='Plot data')


    parser.add_argument('--input', dest='input')
    parser.add_argument('--colvars', dest='colvars',nargs='+')
    args = parser.parse_args()


    sigma=get_sigma(args.input,args.colvars) 
    print(",".join(sigma))


if __name__=="__main__":
    main()