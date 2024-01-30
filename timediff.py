from datetime import datetime
import argparse





def time_diff(start,end):
    # start time
    start_time = start
    end_time = end

    # convert time string to datetime
    t1 = datetime.strptime(start_time, "%Y-%m-%d %H:%M:%S")
    t2 = datetime.strptime(end_time, "%Y-%m-%d %H:%M:%S")
    delta = t2 - t1
    return(t1,t2,delta)

def main():
    parser = argparse.ArgumentParser(
                    prog='TimeDiff',
                    description='Calculate time difference in seconds',
                    epilog='Time is nothing, timing is everything')
    parser.add_argument('--start', type=str, help='Start time (YYYY-MM-DD HH:MM:SS)')
    parser.add_argument('--end', type=str, help='End time (YYYY-MM-DD HH:MM:SS)')
    args = parser.parse_args()
    t1,t2,delta=time_diff(args.start,args.end)
    print('Start time:', t1)
    print('End time:', t2)
    print(f"Time difference is {delta.total_seconds()} seconds")

if __name__=="__main__":
    main()
    
    
