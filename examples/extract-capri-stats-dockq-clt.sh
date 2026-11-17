#!/bin/bash
#

for i in `ls -d $1/*caprieval/capri_clt.tsv`;
do
    export nummodels=`grep -v "#" $i |grep -v clust | wc -l |awk '{print $1}'`
    export naccept=`grep -v "#" $i  |grep -v clust| awk '$13 >= 0.23' | wc -l |awk '{print $1}'`
    export nmedium=`grep -v "#" $i  |grep -v clust| awk '$13 >= 0.49' | wc -l |awk '{print $1}'`
    export nhigh=`grep -v "#" $i  |grep -v clust| awk '$13 >= 0.8' | wc -l |awk '{print $1}'`
    echo "=============================================="
    echo "==" $i
    echo "=============================================="
    echo "Total number of acceptable or better clusters: " $naccept " out of " $nummodels
    echo "Total number of medium or better clusters:     " $nmedium " out of " $nummodels
    echo "Total number of high quality clusters:         " $nhigh " out of " $nummodels
    echo " "
    export firstacceptrank=`grep -v "#" $i |grep -v clust | awk '$13 >= 0.23' | head -1 | awk '{print $1}'`
    export firstacceptirmsd=`grep -v "#" $i |grep -v clust | awk '$13 >= 0.23' | head -1 | awk '{print $7}'`
    export firstacceptfnat=`grep -v "#" $i |grep -v clust | awk '$13 >= 0.23' | head -1 | awk '{print $9}'`
    export firstacceptdockq=`grep -v "#" $i |grep -v clust | awk '$13 >= 0.23' | head -1 | awk '{print $13}'`
    export firstmediumrank=`grep -v "#" $i |grep -v clust | awk '$13 >= 0.49' | head -1 | awk '{print $1}'`
    export firstmediumirmsd=`grep -v "#" $i |grep -v clust | awk '$13 >= 0.49' | head -1 | awk '{print $7}'`
    export firstmediumfnat=`grep -v "#" $i |grep -v clust | awk '$13 >= 0.49' | head -1 | awk '{print $9}'`
    export firstmediumdockq=`grep -v "#" $i |grep -v clust | awk '$13 >= 0.49' | head -1 | awk '{print $13}'`
    export bestrank=`grep -v "#" $i |grep -v clust | sort -nk7 | head -1 | awk '{print $1}'`
    export bestirmsd=`grep -v "#" $i |grep -v clust | sort -nk7 | head -1 | awk '{print $7}'`
    export bestfnat=`grep -v "#" $i |grep -v clust | sort -nk7 | head -1 | awk '{print $9}'`
    export bestdockq=`grep -v "#" $i |grep -v clust | sort -nk7 | head -1 | awk '{print $13}'`
    echo "First acceptable cluster - rank: " $firstacceptrank " i-RMSD: " $firstacceptirmsd " Fnat: " $firstacceptfnat " DockQ: " $firstacceptdockq
    echo "First medium cluster     - rank: " $firstmediumrank " i-RMSD: " $firstmediumirmsd " Fnat: " $firstmediumfnat " DockQ: " $firstmediumdockq
    echo "Best cluster             - rank: " $bestrank " i-RMSD: " $bestirmsd " Fnat: " $bestfnat " DockQ: " $bestdockq 
done
