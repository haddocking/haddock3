#!/bin/csh
#
setenv WDIR /data/capri/Capri47/Target160-scoring/ana_scripts

goto continue

echo ' running analysis...'
$WDIR/run_anal.csh
#
echo ' calculating rmsd matrix and clustering...'
#
continue:

$WDIR/cluster-rmsd.csh
cd rmsds
$HADDOCKTOOLS/cluster_struc target41.rmsd_matrix 5.0 1 >cluster5.0-1.out
$HADDOCKTOOLS/cluster_struc target41.rmsd_matrix 5.0 2 >cluster5.0-2.out
$HADDOCKTOOLS/cluster_struc target41.rmsd_matrix 5.0 4 >cluster5.0-4.out
cd ..

#
echo ' analysing and ranking clusters... min size=4'
#
$WDIR/ana_clusters.csh -best 4 rmsds/cluster5.0-4.out
$WDIR/cluster-scores.csh 4
#
echo ' creating a new cluster directory and copying corresponding files...'
#
mkdir new-clusters-min4
cd new-clusters-min4
$WDIR/copy_selected.csh target41 30 ../clusters_best4.scores 4
#
echo ' running analysis of new clusters...'
#
$WDIR/run_analysis-selected.csh
cd ../
#

#
echo ' analysing and ranking clusters... min size=2'
#
$WDIR/ana_clusters.csh -best 2 rmsds/cluster5.0-2.out
$WDIR/cluster-scores.csh 2
#
echo ' creating a new cluster directory and copying corresponding files...'
#
mkdir new-clusters-min2
cd new-clusters-min2
$WDIR/copy_selected.csh target41 30 ../clusters_best2.scores 2
#
echo ' running analysis of new clusters...'
#
$WDIR/run_analysis-selected.csh
cd ../../
echo ' DONE'

exit:
