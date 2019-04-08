#!/bin/csh

#performs rigid-body pairwise ligand RMSD fitting

#To be used for docking ensembles of two proteins, that have been obtained by rigid body docking (it0)

#This method is not suitable for interface RMSD fitting, but with adaptation it could be applied to
# ligand-interface RMSD (fitting on the receptor, RMSD calculation over the interface residues of the ligand)

#supply:
#- a file containing a list of all structures,
#- a chain to fit on (receptor chain)
#- a chain to calculate the RMSD over (ligand chain)

#There are two steps FILTER1 and FILTER2, where you select atoms. In FILTER1, the atoms used for fitting are
# selected, whereas in FILTER2 the atoms used for RMSD calculation are selected.
# Hence, FILTER2 is a subset of FILTER1
# Currently, FILTER1 is set to all-atom and FILTER2 is set to backbone

#RESULT: pairwise RMSD matrix, printed on screen

if ($# < 3) then
  echo "  Missing arguments"
  echo "  SYNTAX: rbrmsd.csh <list file of structures> <receptor chain> <ligand chain>"
  exit
endif

#set rundir = `echo $0 | awk -v d=$0:h"/" '{i = index($0,"/") } {if (i == 0) $0 = "./"; else $0 = d } {print $0}'`
set rundir = /home/software/haddock/haddock_scripts/smoothing

set d = `echo $1 | awk -v d=$1:h"/" '{i = index($0,"/") } {if (i == 0) $0 = "./"; else $0 = d } {print $0}'`

awk 'NF > 0' $1 > tmp

echo "Copying and filtering structures..." > /dev/stderr
# FILTER1: select all atoms used for fitting
foreach i(`cat tmp`)
  ### CHOOSE ONE OF THESE:###
  
  ### 1: CA filter ###
  pdb_xsegchain $d/$i | awk 'substr($0,27,1) == " "' | awk '$1 == "ATOM" && $3 == "CA" ' | awk '{$0 = substr($0,0,56) "1.00" substr($0,60)} {print $0}' >$i:r.pdb0 
  
  ### 2: backbone filter ###
  #pdb_xsegchain $d/$i | awk 'substr($0,27,1) == " "' | awk '$1 == "ATOM" && ($3 == "CA" || $3 == "C" || $3 == "O" || $3 == "N")' | awk '{$0 = substr($0,0,56) "1.00" substr($0,60)} {print $0}' >$i:r.pdb0 
  
  ### 3: all atom filter ###
  #pdb_xsegchain $d/$i | awk 'substr($0,27,1) == " "' | awk '$1 == "ATOM"' | awk '{$0 = substr($0,0,56) "1.00" substr($0,60)} {print $0}' >$i:r.pdb0 

end

set refe = `head -1 tmp`
set refe = $refe:r.pdb0

echo "Performing initial ProFit fitting..." > /dev/stderr
foreach i(`cat tmp`)  
  set mobi = $i:r.pdb0
    profit $refe $mobi << _Eod_ > out
    zone $2*
    fit
    write $i:r.pdb1
    quit
_Eod_
end

foreach i(`cat tmp`)  
 rm -f $i:r.pdb0
end

echo "Filtering for ligand atoms to be used in RMSD calculation" > /dev/stderr
# FILTER2: select all atoms used for RMSD calculation
foreach i(`cat tmp`)
   ### CHOOSE ONE OF THESE:###
   
   ### 1: CA filter ###
    awk -v c=$3 '$3 == "CA" && substr($0,22,1) == c' $i:r.pdb1 > $i:r.pdb2
  
  ### 2: backbone filter ###
    #awk -v c=$3 '($3 == "CA" || $3 == "C" || $3 == "O" || $3 == "N") && substr($0,22,1) == c' $i:r.pdb1 > $i:r.pdb2
  
  ### 3: all-atom filter ###
    #awk -v c=$3 'substr($0,22,1) == c' $i:r.pdb1 > $i:r.pdb2
  
  rm -f $i:r.pdb1
end

echo "Calculating rotation matrices using ProFit..." > /dev/stderr
set refe = $refe:r.pdb2
cat tmp | wc -l > all.mat
foreach i(`cat tmp`)
  profit $refe $i:r.pdb2 << _Eod_ | awk '$1 == "Mobile" && $2 == "CofG..." {getline; print $0; getline; getline; print $0; getline; print $0; getline; print $0}' >> all.mat
    fit
    matrix
    quit
_Eod_
end

echo "Calculating pairwise RMSD..." > /dev/stderr
$rundir/rbrmsd $refe all.mat

echo "Cleaning up..." > /dev/stderr
foreach i(`cat tmp`)
  rm -f $i:r.pdb2
end

rm -f out
rm -f all.mat
rm -f tmp
