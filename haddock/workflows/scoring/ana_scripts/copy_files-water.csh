#!/bin/csh
#
setenv RDIR /data/capri/Capri47/Target160-scoring

if ($#argv < 3) goto usage

set run=$1
set fileroot=$2
set nbest=$3
set ifile=$4

\cp $RDIR/run$run/run.cns .
\cp $RDIR/run$run/data/new.html .
\cp $RDIR/run$run/data/distances/*.tbl .
\cp $RDIR/run$run/structures/it1/water/structures_haddock-sorted.stat .
\cp $RDIR/run$run/structures/it1/water/structures.scores .
\cp $RDIR/run$run/structures/it1/water/structures.ranking .
\cp $RDIR/run$run/structures/it1/water/structures.ranking-unsorted .
\cp $RDIR/run$run/structures/it1/water/clusters_haddock-sorted.stat_best* .
\cp $RDIR/run$run/structures/it1/water/clusters_nbw-sorted.stat_best* .
\cp $4 .

foreach iclu (`grep -v "#" $ifile |awk '{print $1}' |head -$nbest`)

 set inum=1
   if (-e $RDIR'/run'$run'/structures/it1/water/file.nam_'$iclu'_best4ranked') then
     set filelist=`echo $RDIR'/run'$run'/structures/it1/water/file.nam_'$iclu'_best4ranked'`
     set nc=4
   else if (-e $RDIR'/run'$run'/structures/it1/water/file.nam_'$iclu'_best2ranked') then
     set filelist=`echo $RDIR'/run'$run'/structures/it1/water/file.nam_'$iclu'_best2ranked'`
     set nc=2
   endif
   foreach i ( `head -$nbest $filelist` )
     set oname=`echo $fileroot'_run'$run'_'$iclu'_'$inum'.pdb'`
     grep REMARK $RDIR/run$run/structures/it1/water/$i >$oname
     $HADDOCKTOOLS/pdb_segid-to-chain $RDIR/run$run/structures/it1/water/$i |pdb_selchain -A | grep ATOM >> $oname
     $HADDOCKTOOLS/pdb_segid-to-chain $RDIR/run$run/structures/it1/water/$i |pdb_selchain -C | pdb_chain -A |pdb_xseg -A | grep ATOM >> $oname
     $HADDOCKTOOLS/pdb_segid-to-chain $RDIR/run$run/structures/it1/water/$i |pdb_selchain -B | grep ATOM >> $oname
     $HADDOCKTOOLS/pdb_segid-to-chain $RDIR/run$run/structures/it1/water/$i |pdb_selchain -D | pdb_chain -B |pdb_xseg -B  grep ATOM >> $oname
     echo TER >>$oname
     pdb_selxseg -WAT1 $RDIR/run$run/structures/it1/water/`echo $i|sed s/w/_h2o/` |pdb_chain -C |pdb_xseg -C | pdb_reres -1000 | grep ATOM | sed -e 's/TIP3/WAT\ /' >> $oname
     pdb_selxseg -WAT2 $RDIR/run$run/structures/it1/water/`echo $i|sed s/w/_h2o/` |pdb_chain -C |pdb_xseg -C | pdb_reres -2000 | grep ATOM | sed -e 's/TIP3/WAT\ /'>> $oname
     echo END >>$oname
     @ inum+=1
   end
   \cp $RDIR'/run'$run'/structures/it1/water/file.nam_'$iclu'_best'$nc'ranked' .
   head -1 structures.scores >'file.nam_'$iclu'_best'$nc'ranked-scores'
   foreach i (`cat 'file.nam_'$iclu'_best'$nc'ranked'`)
     grep $i structures.scores >>'file.nam_'$iclu'_best'$nc'ranked-scores'
   end
end

goto exit

usage:

echo 'copy_files.csh: copy best structures of each cluster from an HADDOCK run'
echo ' '
echo ' usage: copy_files.csh run-number fileroot number-of-clusters-to-copy cluster-ranking-filename'
echo ' '
echo ' Note that this script assume that the run directory is  located  in ../../'

exit:
