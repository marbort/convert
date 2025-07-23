sigma=0.025
MIN=$2
MAX=$3
ZMAX=$4
CV=$5
SKIP=1000
stride=$(python ~/SHARED/GitHub/mio/convert/get_stride_colvar.py)
echo "Sigma= $sigma"
head -n 1 colvar
bias=$(python ~/SHARED/GitHub/mio/convert/find_bias.py)

python ~/Downloads/others/opes-metad/FES_from_Reweighting.py --colvar colvar --outfile fes-rew.dat --sigma $sigma --temp 300 --blocks 5 --cv $1 --bias $bias --min $MIN --max $MAX --skiprows $SKIP
python ~/Downloads/others/opes-metad/FES_from_Reweighting.py --colvar colvar --outfile fes-rew.dat --sigma $sigma --temp 300 --stride $stride --cv $1 --bias $bias --min $MIN --max $MAX --skiprows $SKIP
python ~/SHARED/GitHub/mio/convert/plot_single_fep.py --input fes-rew.dat --labx "$CV" --zmax $ZMAX 
