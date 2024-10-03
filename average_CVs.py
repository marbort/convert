import numpy as np

CV1=np.loadtxt('fes-stat_square_sparse.dat_cv1.dat',unpack=True)
CV2=np.loadtxt('fes-stat_square_sparse.dat_cv2.dat',unpack=True)

avg=(CV1[1]+CV2[1])/2
with open('CV_avg.dat','w') as ofile:
   for i,item in enumerate(avg):
        ofile.write(f"{CV1[0][i]} {item}\n")
