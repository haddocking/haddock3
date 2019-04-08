#!/bin/csh
#
setenv RDIR /data/capri/Capri47/Target160-scoring
if ($#argv < 4) goto usage

set run=$1
set fileroot=$2
set nbest=$3
set iclu=$4


set inum=1
set filelist=`echo '../../run'$run'/structures/it1/water/file.nam_clust'$iclu`
set oclus=`echo $fileroot'_run'$run'_clust'$iclu'.scores'`
\cp $RDIR'/run'$run'/structures/it1/water/file.nam_clust'$iclu'_best'$nbest'-scores' $oclus
foreach i ( `head -$nbest $filelist` )
     set oname=`echo $fileroot'_run'$run'_clust'$iclu'_'$inum'.pdb'`
     grep REMARK $RDIR/run$run/structures/it1/water/$i >$oname
     $HADDOCKTOOLS/pdb_segid-to-chain $RDIR/run$run/structures/it1/water/$i |pdb_selchain -A | grep ATOM >> $oname
     $HADDOCKTOOLS/pdb_segid-to-chain $RDIR/run$run/structures/it1/water/$i |pdb_selchain -B | grep ATOM >> $oname
     $HADDOCKTOOLS/pdb_segid-to-chain $RDIR/run$run/structures/it1/water/$i |pdb_selchain -C | pdb_chain -A | grep ATOM >> $oname
     $HADDOCKTOOLS/pdb_segid-to-chain $RDIR/run$run/structures/it1/water/$i |pdb_selchain -D | pdb_chain -B | grep ATOM >> $oname
     echo END >>$oname
     @ inum+=1
end


goto exit

usage:

echo 'copy_files.csh: copy best structures of a specific cluster from an HADDOCK run'
echo ' '
echo ' usage: copy_files.csh run-number fileroot number-of-structures-to-copy cluster-number'
echo ' '
echo ' Edit in this script the location of the run directory'

exit:
