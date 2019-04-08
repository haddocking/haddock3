#!/bin/csh
#
if ($#argv < 3) goto usage

set fileroot=$1
set nbest=$2
set ifile=$3
set ibest=$4

foreach iclu (`grep -v "#" $ifile |awk '{print $1}' |head -$nbest`)

   set inum=1

#   set filelist=`echo '../file.nam_'$iclu'_best'$ibest'ranked'`
   set filelist=`echo '../file.nam_'$iclu'_best'$ibest`
   cp $filelist'-scores' .
   foreach i ( `head -$ibest $filelist` )
     set oname=`echo $fileroot'_'$iclu'_'$inum'.pdb'`
     sed -e 's/TIP3B/TIP3C/' ../$i > $oname
     @ inum+=1
   end

end
\cp ../structures_haddock-sorted.stat .
\cp ../structures.scores .
\cp ../structures.ranking .
\cp ../structures.ranking-unsorted .
\cp ../clusters_haddock-sorted.stat_best* .
\cp ../clusters_nbw-sorted.stat_best* .
\cp $3 .

goto exit

usage:

echo 'copy_files.csh: copy best structures of each cluster from an HADDOCK run'
echo ' '
echo ' usage: copy_files.csh run-number fileroot number-of-clusters-to-copy cluster-ranking-filename nbest'
echo ' '
echo ' Note that this script assume that the files to copy are located in ../'

exit:
