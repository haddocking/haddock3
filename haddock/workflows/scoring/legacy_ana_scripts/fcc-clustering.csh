#!/bin/csh
#
setenv WDIR /home/joao/Code/fcc

python  $WDIR/scripts/make_contacts.py -n4 -e $WDIR/src/contact_fcc -f $1

sed -e 's/pdb/contacts/' $1 |grep contacts > file.contacts

python $WDIR/scripts/calc_fcc_matrix.py -f file.contacts -o fcc_matrix.out

python $WDIR/scripts/cluster_fcc.py fcc_matrix.out 0.75 -o cluster_0.75-4.out -c 4
python $WDIR/scripts/cluster_fcc.py fcc_matrix.out 0.75 -o cluster_0.75-2.out -c 2


