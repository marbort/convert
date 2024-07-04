#! /usr/bin/zsh

FILE=$1
TOPO=$2
AROUND=$3

echo " parm $TOPO \n trajin $FILE ucell lastframe \n trajout ${FILE}.pdb \n go \n quit " > cpptraj.script

cpptraj -i cpptraj.script

python ~/SHARED/GitHub/mio/convert/last_dcd_to_pdb.py $TOPO ${FILE}.pdb $AROUND



