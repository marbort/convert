from datetime import datetime
import argparse

parser = argparse.ArgumentParser(
                    prog='TimeDiff',
                    description='Calculate time difference in seconds',
                    epilog='Time is nothing, timing is everything')

  
parser.add_argument('--start', type=str, help='Start time (HH:MM:SS)')
parser.add_argument('--end', type=str, help='End time (HH:MM:SS)')

args = parser.parse_args()


# start time
start_time = args.start
end_time = args.end

# convert time string to datetime
t1 = datetime.strptime(start_time, "%H:%M:%S")
print('Start time:', t1.time())

t2 = datetime.strptime(end_time, "%H:%M:%S")
print('End time:', t2.time())

# get difference
delta = t2 - t1

if delta.total_seconds() < 0:
    print(f"Time difference is {86400 + delta.total_seconds()} seconds")

else:
    # time difference in seconds
    print(f"Time difference is {delta.total_seconds()} seconds")

