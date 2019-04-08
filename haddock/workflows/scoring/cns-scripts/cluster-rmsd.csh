#!/bin/csh
#

rm -rf rmsds/
mkdir rmsds

set inum=1

foreach i (`cat file.nam`)
  $HADDOCKTOOLS/pdb_segid-to-chain $i |grep ' CA '|zone_select -zA167,A430 >rmsds/target26_$inum.pdb
  $HADDOCKTOOLS/pdb_segid-to-chain $i |grep ' CA '|zone_select -zB45,B151 >>rmsds/target26_$inum.pdb
  echo END >>rmsds/target26_$inum.pdb
  @ inum+=1
end

set numstr=`wc -l file.nam |awk '{print $1}'`
cd rmsds
$HADDOCKTOOLS/make_rmsd_matrix.csh A B $numstr target26

