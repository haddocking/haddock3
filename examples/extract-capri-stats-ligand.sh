#!/bin/bash
#

for i in `ls -d $1/*caprieval/capri_ss.tsv`;
do
    export nummodels=`grep pdb $i | wc -l |awk '{print $1}'`
    export naccept=`grep pdb $i | awk '$8 <= 2.5' | wc -l |awk '{print $1}'`
    export nmedium=`grep pdb $i | awk '$8 <= 1' | wc -l |awk '{print $1}'`
    export nhigh=`grep pdb $i | awk '$8 <= 0.5' | wc -l |awk '{print $1}'`
    echo "=============================================="
    echo "== Assessing based in ilRMSD with ligand criteria"
    echo "==" $i
    echo "=============================================="
    echo "Total number of acceptable or better models: " $naccept " out of " $nummodels
    echo "Total number of medium or better models:     " $nmedium " out of " $nummodels
    echo "Total number of high quality models:         " $nhigh " out of " $nummodels
    echo " "
    export firstacceptrank=`grep pdb $i | awk '$8 <= 2.5' | head -1 | awk '{print $3}'`
    export firstacceptirmsd=`grep pdb $i | awk '$8 <= 2.5' | head -1 | awk '{print $8}'`
    export firstacceptfnat=`grep pdb $i | awk '$8 <= 2.5' | head -1 | awk '{print $6}'`
    export firstacceptdockq=`grep pdb $i | awk '$8 <= 2.5' | head -1 | awk '{print $9}'`
    export firstmediumrank=`grep pdb $i | awk '$8 <= 1.0' | head -1 | awk '{print $3}'`
    export firstmediumirmsd=`grep pdb $i | awk '$8 <= 1.0' | head -1 | awk '{print $8}'`
    export firstmediumfnat=`grep pdb $i | awk '$8 <= 1.0' | head -1 | awk '{print $6}'`
    export firstmediumdockq=`grep pdb $i | awk '$8 <= 1.0' | head -1 | awk '{print $9}'`
    export bestrank=`grep pdb $i | sort -nk8 | head -1 | awk '{print $3}'`
    export bestirmsd=`grep pdb $i | sort -nk8 | head -1 | awk '{print $8}'`
    export bestfnat=`grep pdb $i | sort -nk8 | head -1 | awk '{print $6}'`
    export bestdockq=`grep pdb $i | sort -nk8 | head -1 | awk '{print $9}'`
    echo "First acceptable model - rank: " $firstacceptrank " il-RMSD: " $firstacceptirmsd " Fnat: " $firstacceptfnat " DockQ: " $firstacceptdockq
    echo "First medium model     - rank: " $firstmediumrank " il-RMSD: " $firstmediumirmsd " Fnat: " $firstmediumfnat " DockQ: " $firstmediumdockq
    echo "Best model             - rank: " $bestrank " il-RMSD: " $bestirmsd " Fnat: " $bestfnat " DockQ: " $bestdockq
done
