sigma=$(tail -n1 KERNELS | awk '{printf "%.3f", $3}')
echo "Sigma= $sigma"
head -n 1 colvar
echo "Select bias(es) to use (comma separated)"
read bias

python ~/Downloads/others/opes-metad/FES_from_Reweighting.py --colvar colvar --outfile fes-rew.dat --sigma $sigma --temp 300 --blocks 5 --cv $1 --bias $bias
