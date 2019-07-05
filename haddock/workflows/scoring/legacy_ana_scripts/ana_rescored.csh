#!/bin/csh
#
source /home/abonvin/haddock2.4/haddock_configure.csh
setenv WDIR /data/capri/Capri47/Target160-scoring/ana_scripts

goto continue

foreach i (*conv.pdb)
   $HADDOCKTOOLS/pdb_segid-to-chain $i >ttt
   \mv ttt $i
end

echo ' Sorting and generating file.nam file.list...'
$WDIR/make-files-init.csh *[01]_conv.pdb
$WDIR/make-files-cont.csh *[23]_conv.pdb
$WDIR/make-files-cont.csh *[45]_conv.pdb
$WDIR/make-files-cont.csh *[67]_conv.pdb
$WDIR/make-files-end.csh *[89]_conv.pdb
#
echo ' Now running analysis of single structures...'
#$WDIR/run_anal-it0.csh
$WDIR/run_anal.csh
#
#
echo ' calculating fcc matrix and clustering...'
#
mkdir rmsds
python $HADDOCKTOOLS/make_contacts.py `cat file.nam` --exec $HADDOCKTOOLS/contact_fcc
cat file.nam |sed -e 's/pdb/contacts/' >file.contacts
python $HADDOCKTOOLS/calc_fcc_matrix.py -f file.contacts -o rmsds/target160.fcc_matrix -i >& rmsds/fcc.out

cd rmsds
python $HADDOCKTOOLS/cluster_fcc.py target160.fcc_matrix 0.75 -c 1  >cluster0.75-1.out
python $HADDOCKTOOLS/cluster_fcc.py target160.fcc_matrix 0.75 -c 2  >cluster0.75-2.out
python $HADDOCKTOOLS/cluster_fcc.py target160.fcc_matrix 0.75 -c 4  >cluster0.75-4.out
python $HADDOCKTOOLS/cluster_fcc.py target160.fcc_matrix 0.6  -c 1  >cluster0.6-1.out
python $HADDOCKTOOLS/cluster_fcc.py target160.fcc_matrix 0.6  -c 2  >cluster0.6-2.out
python $HADDOCKTOOLS/cluster_fcc.py target160.fcc_matrix 0.6  -c 4  >cluster0.6-4.out
cd ..

#
echo ' analysing and ranking clusters... min size=4'
#
$WDIR/ana_clusters.csh -best 4 rmsds/cluster0.6-4.out >&/dev/null

$WDIR/cluster-scores.csh 4
#
echo ' creating a new cluster directory and copying corresponding files...'
#
mkdir new-clusters-min4
cd new-clusters-min4
$WDIR/copy_selected.csh target160 30 ../clusters_best4.scores 4
#
echo ' running analysis of new clusters...'
#
$WDIR/run_analysis-selected.csh
cd ../
#
#goto exit
#
echo ' analysing and ranking clusters... min size=2'
#
$WDIR/ana_clusters.csh -best 2 rmsds/cluster0.6-2.out >&/dev/null
$WDIR/cluster-scores.csh 2
#
echo ' creating a new cluster directory and copying corresponding files...'
#
mkdir new-clusters-min2
cd new-clusters-min2
$WDIR/copy_selected.csh target160 30 ../clusters_best2.scores 2
#
echo ' running analysis of new clusters...'
#
$WDIR/run_analysis-selected.csh
cd ../../

continue:

echo ' analysing and ranking clusters... min size=1'
#
$WDIR/ana_clusters.csh -best 2 rmsds/cluster0.6-1.out >&/dev/null
$WDIR/cluster-scores.csh 1
#
echo ' creating a new cluster directory and copying corresponding files...'
#
mkdir new-clusters-min1
cd new-clusters-min1
$WDIR/copy_selected.csh target160 30 ../clusters_best1.scores 1
#
echo ' running analysis of new clusters...'
#
$WDIR/run_analysis-selected.csh
cd ../../

exit:
echo ' DONE'
