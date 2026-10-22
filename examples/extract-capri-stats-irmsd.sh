#!/bin/bash
#

for i in `ls -d $1/*caprieval/capri_ss.tsv`;
do
    export nummodels=`grep pdb $i | wc -l |awk '{print $1}'`
    export naccept=`grep pdb $i | awk '$5 <= 4' | wc -l |awk '{print $1}'`
    export nmedium=`grep pdb $i | awk '$5 <= 2' | wc -l |awk '{print $1}'`
    export nhigh=`grep pdb $i | awk '$5 <= 1' | wc -l |awk '{print $1}'`
    echo "=============================================="
    echo "==" $i
    echo "=============================================="
    echo "Total number of acceptable or better models: " $naccept " out of " $nummodels
    echo "Total number of medium or better models:     " $nmedium " out of " $nummodels
    echo "Total number of high quality models:         " $nhigh " out of " $nummodels
    echo " "
    export firstacceptrank=`grep pdb $i | awk '$5 <= 4' | head -1 | awk '{print $3}'`
    export firstacceptirmsd=`grep pdb $i | awk '$5 <= 4' | head -1 | awk '{print $5}'`
    export firstacceptfnat=`grep pdb $i | awk '$5 <= 4' | head -1 | awk '{print $6}'`
    export firstacceptdockq=`grep pdb $i | awk '$5 <= 4' | head -1 | awk '{print $9}'`
    export firstmediumrank=`grep pdb $i | awk '$5 <= 2' | head -1 | awk '{print $3}'`
    export firstmediumirmsd=`grep pdb $i | awk '$5 <= 2' | head -1 | awk '{print $5}'`
    export firstmediumfnat=`grep pdb $i | awk '$5 <= 2' | head -1 | awk '{print $6}'`
    export firstmediumdockq=`grep pdb $i | awk '$5 <= 2' | head -1 | awk '{print $9}'`
    export bestrank=`grep pdb $i | sort -nk5 | head -1 | awk '{print $3}'`
    export bestirmsd=`grep pdb $i | sort -nk5 | head -1 | awk '{print $5}'`
    export bestfnat=`grep pdb $i | sort -nk5 | head -1 | awk '{print $6}'`
    export bestdockq=`grep pdb $i | sort -nk5 | head -1 | awk '{print $9}'`
    echo "First acceptable model - rank: " $firstacceptrank " i-RMSD: " $firstacceptirmsd " Fnat: " $firstacceptfnat " DockQ: " $firstacceptdockq
    echo "First medium model     - rank: " $firstmediumrank " i-RMSD: " $firstmediumirmsd " Fnat: " $firstmediumfnat " DockQ: " $firstmediumdockq
    echo "Best model             - rank: " $bestrank " i-RMSD: " $bestirmsd " Fnat: " $bestfnat " DockQ: " $bestdockq
done
