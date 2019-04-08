#!/bin/csh
#
if ($#argv < 3) goto usage
#
# Define the location of profit
#
if ( `printenv |grep PROFIT | wc -l` == 0) then
  set found=`which profit |wc -l`
  if ($found == 0) then
     echo 'PROFIT environment variable not defined'
     echo 'and not profit not found on your system'
     echo ' ==> stopping'
     goto exit
  else
     setenv PROFIT `which profit`
  endif
endif

set i=0
set inum=$3
set jnum=$3
@ inum-=1

cat /dev/null >rmsd_matrix.profit
while ($i < $inum) 
  @ i+=1
  set iname=$4_$i.pdb
  echo refe $iname >>rmsd_matrix.profit
  set j=$i
  while ($j < $jnum)
    @ j+=1
    set jname=$4_$j.pdb
    echo 'zone B100-B295' >>rmsd_matrix.profit
    echo atom CA,C,N  >>rmsd_matrix.profit
    echo mobi $jname >>rmsd_matrix.profit
    echo fit >>rmsd_matrix.profit
    echo rzone A77-A180 >>rmsd_matrix.profit
    echo rzone A200-A270 >>rmsd_matrix.profit
    echo zone clear >>rmsd_matrix.profit
  end
end
echo quit >>rmsd_matrix.profit

$PROFIT <rmsd_matrix.profit>&rmsd_matrix.out

grep RMS rmsd_matrix.out | awk 'BEGIN{count=3} {count=count-1; if (! count) {print $2;count=3}}' >rmsd_matrix.tmp

cat rmsd_matrix.tmp | awk -v totstruc=$3 'BEGIN{i=1;j=2} {print i,j,$1;j=j+1;if(j>totstruc) {i=i+1;j=i+1}}'>$4.rmsd_matrix
#\rm rmsd_matrix.out rmsd_matrix.tmp

goto exit

usage:

echo 'make_rmsd_matrix.csh:  calculates the pairwise backbone RMSD matrix for clustering'
echo '                       fitting on defined chain, RMSD calculation on complete backbone'
echo ' '
echo '    usage: make_rmsd_matrix.csh chainID-for-fitting chainID-for-rmsd Number_of_structures filename_root'
echo ' '
echo '    The chainID defines which chain will be used for  superposition'
echo '    The PDB filenames are assumed to be: filename_root_#.pdb '
echo '    where # indicates the structure number'

exit:
