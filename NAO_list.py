import sys
seq=[]
for i in range(int(sys.argv[1]),int(sys.argv[2])+1):
    seq.append(str(i))
print(seq)
print(",".join(seq))