import numpy as np
import argparse





def extract_data(input,cv1_val,cv2_val,interval,TS=False):
    points=[]
    min_val=10000
    
    max_val=0
    
    data=np.loadtxt(input)
    for i in data:
        if  i[0]-float(interval) < float(cv1_val) < i[0]+float(interval) and i[1]-float(interval)< float(cv2_val) < i[1]+float(interval):
            points.append(i)
    
    for i in points:
        if i[-1] < min_val:
            min_val=i[-1]
            min_arr=i
        if i[-1] > max_val:
            max_val=i[-1]
            max_arr=i
            

    if TS:
        print(f"TS at {max_arr[0]}, {max_arr[1]} energy: {max_arr[2]}") 
    else:
        print(f"Min at {min_arr[0]}, {min_arr[1]} energy: {min_arr[2]}")
    return(min_val)

def main(): 
    parser = argparse.ArgumentParser(description='Extract energy from FES at given CV values')
    parser.add_argument('--input', type=str, help='Input file',default=None)
    parser.add_argument('--cv1', type=float, help='Value of CV1',default=0.0)
    parser.add_argument('--cv2', type=float, help='Value of CV2',default=0.0)
    parser.add_argument('--interval', type=float, help='Interval around CV values to consider',default=0.2)
    parser.add_argument('--TS', action='store_true', help='Extract transition state (max) instead of min',default=False)
    
    args = parser.parse_args()
    
    extract_data(args.input,args.cv1,args.cv2,args.interval,args.TS)
if __name__=="__main__":
    main()