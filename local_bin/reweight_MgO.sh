#!/bin/sh

CV1=$(grep opes plumed.dat | awk '{print $3}' | cut -d "=" -f2 | cut -d "," -f1)
CV2=$(grep opes plumed.dat | awk '{print $3}' | cut -d "=" -f2 | cut -d "," -f2)
min1=-0.5
max1=4.5
min2=-0.5
max2=4.5
sigma1=$(tail -n1 KERNELS | awk '{printf "%.3f", $4}')
sigma2=$(tail -n1 KERNELS | awk '{printf "%.3f", $5}')
temp=300
skip=50
stride=$(grep " dump_freq " start.lmp | cut -d " " -f11)
driver_cmd="plumed driver --plumed plumed_driver.dat --ixyz traj.xyz --timestep 0.00025 --trajectory-stride $stride --length-units A"
typ_elm="/run/media/marco/T7 Shield/SHARED/RATIO/WP1/ML/fox/MeMgCl/UMBRELLA/Mg-O_coord/types_to_elms.json"
file=colvar
minima='structure_minima.dat'
#awk 'NR % 50 == 1' colvar > $file 

echo "Sigma 1 = $sigma1"
echo "Sigma 2 = $sigma2"

params=${CV1}_${CV2}_${min1}_${min2}_${max1}_${max2}_${sigma1}_${sigma2}_${temp}_${skip}.params

if [ -f $params ]
then
 echo "Reweighting done"
else
  python ~/Downloads/others/opes-metad/FES_from_Reweighting.py  --colvar $file --outfile fes-rew_square_sparse_cv1.dat --sigma $sigma1 --temp $temp --cv $CV1 --bias opes.bias --min " $min1" --max $max1 --skiprows $skip
  #
  python ~/Downloads/others/opes-metad/FES_from_Reweighting.py  --colvar $file --outfile fes-rew_square_sparse_cv2.dat --sigma $sigma2 --temp $temp --cv $CV2 --bias opes.bias --min " $min2" --max $max2 --skiprows $skip
  #
  python ~/Downloads/others/opes-metad/FES_from_Reweighting.py  --colvar $file --outfile fes-rew_square_sparse.dat --sigma ${sigma1},$sigma2 --temp $temp --cv $CV1,$CV2 --bias opes.bias --min " ${min1}"," $min2" --max ${max1},$max2 --skiprows $skip

fi

sed "s/opes/#opes/" plumed.dat > plumed_driver.dat
sed -i "s/STRIDE=200/STRIDE=1/" plumed_driver.dat
sed -i "s/FILE=colvar/FILE=colvar_trj/" plumed_driver.dat


##REWEIGHTING###
python ~/SHARED/GitHub/mio/convert/plot2d_reweight.py fes-rew_square_sparse.dat $CV1 $CV2
python ~/SHARED/GitHub/mio/convert/lammpstrj_to_xyzbox.py --input lmp.lammpstrj --output traj.xyz --types_elm "$typ_elm"
#plumed driver --plumed plumed_driver.dat --ixyz traj.xyz --timestep 0.00025 --trajectory-stride "$stride" --length-units A
echo $driver_cmd > plumed_driver.sh
. $PWD/plumed_driver.sh
python ~/SHARED/GitHub/mio/convert/find_CV_struct.py --input colvar_trj --CVs $CV1 $CV2 --val minima.dat --tol 0.1 0.1
python ~/SHARED/GitHub/mio/convert/plot2d_reweight.py fes-rew_square_sparse.dat $CV1 $CV2 $minima
python ~/SHARED/GitHub/mio/convert/plot_1D-2D_seaborn.py fes-rew_square_sparse.dat $CV1 $CV2 fes-rew_square_sparse_cv1.dat fes-rew_square_sparse_cv2.dat $1



rm *.params
touch ${params}
