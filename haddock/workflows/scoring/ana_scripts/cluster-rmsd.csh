#!/bin/csh
#

rm -rf rmsds/
mkdir rmsds

set inum=1

foreach i (`cat file.nam`)
  $HADDOCKTOOLS/pdb_segid-to-chain $i |grep ' CA ' >rmsds/target41_$inum.pdb
  $HADDOCKTOOLS/pdb_segid-to-chain $i |grep ' CA ' >>rmsds/target41_$inum.pdb
  echo END >>rmsds/target41_$inum.pdb
  @ inum+=1
end

set numstr=`wc -l file.nam |awk '{print $1}'`
cd rmsds
../make_rmsd_matrix.csh A B $numstr target41

