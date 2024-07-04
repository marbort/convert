sigma=$(tail -n1 KERNELS | awk '{printf "%.3f", $3}')
echo "Sigma= $sigma"
python ~/Downloads/others/opes-metad/FES_from_Reweighting.py --colvar colvar --outfile fes-rew.dat --sigma $sigma --temp 300 --blocks 5 --cv CV1 --bias lwall.bias,restraint.bias,restraint2.bias,opes.bias --skip 2000
python ~/Downloads/others/opes-metad/FES_from_State.py --temp 300 -f STATE
python ~/SHARED/GitHub/mio/convert/plot_colvars_opes.py 

