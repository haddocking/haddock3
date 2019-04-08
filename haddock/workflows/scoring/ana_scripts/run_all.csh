#!/bin/csh
#
setenv target complex
setenv WDIR /data/capri/Capri47/Target160-scoring/ana_scripts

cd $1/structures/it1/water/analysis
\rm *out.gz; $HADDOCKTOOLS/cluster_struc $target'_rmsd.disp' 7.5 4 >cluster7.5.out
cd ..
$WDIR/run_anal-it0.csh
#$WDIR/fraction-native.csh `cat file.nam`
$WDIR/ana_clusters.csh -best 4 analysis/cluster7.5.out
$WDIR/cluster-scores.csh 4
#$WDIR/cluster-fnat.csh 4
cd ../analysis
\rm *out.gz; $HADDOCKTOOLS/cluster_struc $target'_rmsd.disp' 7.5 4 >cluster7.5.out
cd ..
$WDIR/run_anal-it0.csh
#$WDIR/fraction-native.csh `cat file.nam`
$WDIR/ana_clusters.csh -best 4 analysis/cluster7.5.out
$WDIR/cluster-scores.csh 4
#$WDIR/cluster-fnat.csh 4
cd ../it0
$WDIR/run_anal-it0.csh
set i1=`wc -l file.nam`
@ i1-=5000
#$WDIR/fraction-native.csh `head -5000 file.nam`
#$WDIR/fraction-native-cont.csh `tail -$i1 file.nam`
cd ../../../





