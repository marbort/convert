import numpy as np
import json
import argparse


def find_min(data):
    data=np.array(data)
    min_val=data[0]
    min_index=0
    for i in data[1:]:
        if i<min_val:
            min_val=i
        else:
            break
    min_index=np.where(data==min_val)[0][0]
    print(f"Min value: {min_val} at index {min_index}") 
    return(min_val,min_index)
            
def find_max(data):
    data=np.array(data)
    max_val=data[0]
    max_index=0
    print(max_val)
    for i in data[1:]:
        if i>max_val:
            max_val=i
        else:
            break
    max_index=np.where(data==max_val)[0][0]
    print(f"Max value: {max_val} at index {max_index}") 
    return(max_val,max_index)




def calc_DGs(data,fac_free=1,fac_cv=1):
    with open(data,'r') as data:
        vals=json.load(data)
    free=np.array(vals["summary"]["equil"]["free"])
    cv=np.array(vals["summary"]["equil"]["cv"])
    errors=np.array(vals["summary"]["std_free_mean"])
    prod,prod_idx=find_min(free)
    TS,TS_idx=find_max(free[prod_idx:])
    TS_idx=np.where(free==TS)[0][0] 
    react=np.min(free[TS_idx:])
    react_idx=np.where(free==react)[0][0]
    DGr=(prod-react)*fac_free
    DGact=(TS-react)*fac_free
    cv_prod=cv[prod_idx]
    cv_react=cv[react_idx]
    cv_TS=cv[TS_idx]
    
    react_err=errors[react_idx]
    TS_err=errors[TS_idx]
    prod_err=errors[prod_idx]
    
    print(react_err,TS_err,prod_err)
        
    if react > 0:
        TS=TS-react
        prod=prod-react
        react=0.0
        
    
    return(DGr,DGact,prod*fac_free,TS*fac_free,react*fac_free,cv_prod,cv_react,cv_TS,react_err,prod_err,TS_err)

def calc_DGs_wham(data,fac_free=1,fac_cv=1):
    vals=np.loadtxt(data,unpack=True)
    free=vals[1]*fac_free
    cv=vals[0]*fac_cv
    errors=vals[2]*fac_free
    prod,prod_idx=find_min(free)
    TS,TS_idx=find_max(free[prod_idx:])
    TS_idx=np.where(free==TS)[0][0] 
    react=np.min(free[TS_idx:])
    react_idx=np.where(free==react)[0][0]
    DGr=(prod-react)
    DGact=(TS-react)
    cv_prod=cv[prod_idx]
    cv_react=cv[react_idx]
    cv_TS=cv[TS_idx]
    
    react_err=errors[react_idx]
    TS_err=errors[TS_idx]
    prod_err=errors[prod_idx]
    
    print(react_err,TS_err,prod_err)
        
    if react > 0:
        TS=TS-react
        prod=prod-react
        react=0.0
        
    
    return(DGr,DGact,prod,TS,react,cv_prod,cv_react,cv_TS,react_err,prod_err,TS_err)
    
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("data",type=str,help="json file containing free energy data")
    parser.add_argument("--fac_free",type=float,default=1,help="Factor to multiply the free energy")
    parser.add_argument("--fac_cv",type=float,default=1,help="Factor to multiply the CV")
    parser.add_argument("--multi",action="store_true",help="Read paths from a json file")
    parser.add_argument("--wham",action="store_true",help="Read data from a wham file")
    
    args=parser.parse_args()
    if args.multi:
        results=[]
        with open(args.data,'r') as data:
            a=json.load(data)
        for system in a:
            data=a[system]
            #try:
            if args.wham:
                DGr,DGact,e_prod,e_TS,e_react,cv_prod,cv_react,cv_TS,react_err,prod_err,TS_err=calc_DGs_wham(data,args.fac_free,args.fac_cv)
            else:
                DGr,DGact,e_prod,e_TS,e_react,cv_prod,cv_react,cv_TS,react_err,prod_err,TS_err=calc_DGs(data,args.fac_free,args.fac_cv)
            results.append(f"{system},{cv_react:6.2f},{cv_TS:6.2f},{cv_prod:6.2f},{e_react:6.2f},{react_err:6.2f},{e_TS:6.2f},{TS_err:6.2f},{e_prod:6.2f},{prod_err:6.2f}")
            print(f"{system} OK")
            #except:
            #    print(f"Error in {system}")
        with open("results.csv",'w') as ofile:
            ofile.write("System,REACT_CV,TS_CV,PROD_CV,REACT_ENERGY,REACT_ERR,TS_ENERGY,TS_ERR,PROD_ENERGY,PROD_ERR\n")
            for line in results:
                ofile.write(line)
                ofile.write("\n")
    else:
        data = args.data
        if args.wham:
                DGr,DGact,e_prod,e_TS,e_react,cv_prod,cv_react,cv_TS,react_err,prod_err,TS_err=calc_DGs_wham(data,args.fac_free,args.fac_cv)
    #print format REACT_CV,TS_CV,PROD_CV,REACT_ENERGY,TS_ENERGY,PROD_ENERGY,DGr,DGact
        else:
                DGr,DGact,e_prod,e_TS,e_react,cv_prod,cv_react,cv_TS,react_err,prod_err,TS_err=calc_DGs(data,args.fac_free,args.fac_cv)
        print(f"{cv_react:6.2f},{cv_TS:6.2f},{cv_prod:6.2f},{e_react:6.2f},{react_err:6.2f},{e_TS:6.2f},{TS_err:6.2f},{e_prod:6.2f},{prod_err:6.2f}")
    

if __name__=="__main__":
    main()