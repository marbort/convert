import argparse


def plus(a,b):
    return(a+b)
def minus(a,b):
    return(a-b)


parser = argparse.ArgumentParser(description='Test')
parser.add_argument('--funct', dest='funct', 
                    type=str, help='select function',nargs='+')

args = parser.parse_args()

dict={"plus":plus(1,2),"minus":minus(1,2)}


for i in args.funct:
    a=dict[i]
    print(a)
