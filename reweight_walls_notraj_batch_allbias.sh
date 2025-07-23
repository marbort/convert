#!/bin/bash

CV1=$(grep opes: plumed.dat | awk '{print $3}' | cut -d "=" -f2 | cut -d "," -f1)
CV2=$(grep opes: plumed.dat | awk '{print $3}' | cut -d "=" -f2 | cut -d "," -f2)
min1=-0.1
max1=2.1
min2=-0.1
max2=3.0
#sigma1=$(tail -n1 KERNELS | awk '{printf "%.3f", $4}')
#sigma2=$(tail -n1 KERNELS | awk '{printf "%.3f", $5}')
sigma1=0.025
sigma2=0.025
tot_time=$(tail -n1 colvar | awk '{print $1}')
#sigma1=0.05
#sigma2=0.05
temp=300
skip=1000
stride=$(grep " dump_freq " start.lmp | cut -d " " -f11)
stride_conv=$(python ~/SHARED/GitHub/mio/convert/get_stride_colvar.py)
driver_cmd="plumed driver --plumed plumed_driver.dat --ixyz traj.xyz --timestep 0.00025 --trajectory-stride $stride --length-units A"
typ_elm="/run/media/marco/T7 Shield/SHARED/RATIO/WP1/ML/fox/MeMgCl/UMBRELLA/Mg-O_coord/types_to_elms.json"
file=colvar
minima='structure_minima.dat'
#awk 'NR % 50 == 1' colvar > $file 

params=${tot_time}_${CV1}_${CV2}_${min1}_${min2}_${max1}_${max2}_${sigma1}_${sigma2}_${temp}_${skip}.params
bias=$(python ~/SHARED/GitHub/mio/convert/find_bias.py)

if [ -f $params ]
then
 echo "Reweighting done"
else
  python ~/Downloads/others/opes-metad/FES_from_Reweighting.py  --colvar $file --outfile fes-rew_square_sparse_cv1_walls.dat --sigma $sigma1 --temp $temp --cv $CV1 --bias $bias --min " $min1" --max $max1 --skiprows $skip  
  python ~/Downloads/others/opes-metad/FES_from_Reweighting.py  --colvar $file --outfile fes-rew_square_sparse_cv1_walls.dat --sigma $sigma1 --temp $temp --cv $CV1 --bias $bias --min " $min1" --max $max1 --skiprows $skip --stride $stride_conv 
  #
  python ~/Downloads/others/opes-metad/FES_from_Reweighting.py  --colvar $file --outfile fes-rew_square_sparse_cv2_walls.dat --sigma $sigma2 --temp $temp --cv $CV2 --bias $bias --min " $min2" --max $max2 --skiprows $skip
  #
  python ~/Downloads/others/opes-metad/FES_from_Reweighting.py  --colvar $file --outfile fes-rew_square_sparse_walls.dat --sigma ${sigma1},$sigma2 --temp $temp --cv $CV1,$CV2 --bias $bias --min " ${min1}"," $min2" --max ${max1},$max2 --skiprows $skip
  if [ -f STATE ]
  then
    echo "STATE file exists"
  else
    python ~/Downloads/others/opes-metad/State_from_Kernels.py --kernels KERNELS --outfile STATE
  fi
python ~/Downloads/others/opes-metad/FES_from_State.py  --state STATE --outfile fes-stat_square_sparse_walls.dat --temp $temp  --min " ${min1}"," $min2" --max ${max1},$max2 


fi

python ~/SHARED/GitHub/mio/convert/plot2d_reweight.py --input fes-rew_square_sparse_walls.dat --xlab CV1 --ylab CV2 --max 200 
python ~/SHARED/GitHub/mio/convert/plot_1D-2D_seaborn.py fes-rew_square_sparse_walls.dat CV1 CV2 fes-rew_square_sparse_cv1_walls.dat fes-rew_square_sparse_cv2_walls.dat 200

##STATES##
python ~/SHARED/GitHub/mio/convert/plot2d_reweight.py --input fes-stat_square_sparse_walls.dat --xlab CV1 --ylab CV2  --max 200
python ~/SHARED/GitHub/mio/convert/plot_1D-2D_seaborn.py fes-stat_square_sparse_walls.dat CV1 CV2 fes-stat_square_sparse_walls.dat_cv1.dat fes-stat_square_sparse_walls.dat_cv2.dat 200



rm *.params
touch ${params}
