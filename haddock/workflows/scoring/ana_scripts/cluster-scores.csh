#!/bin/csh
#
if ($#argv == 0) then
  echo "Note that you can specify a number to perform "
  echo "the analysis for the Nbest structures of a cluster"
  echo "... scoring of clusters running..."
else
  echo "... scoring of clusters running including "$1" best structures of each cluster..."
endif

foreach i (file.nam_clust[1-9] file.nam_clust[1-9][0-9] file.nam_clust[1-9][0-9][0-9] file.nam_clust[1-9]_best{$1} file.nam_clust[1-9][0-9]_best{$1} file.nam_clust[1-9][0-9][0-9]_best{$1} )
  head -1 structures.ranking >$i'-scores'
  foreach j (`cat $i`)
    grep $j structures.ranking >>$i'-scores'
  end
end

echo "#cluster size rmsd_emin sd haddock-score sd bsa sd Edesolv sd fc-elec sd fc-desol sd Dfire sd probe sd probe/A**2 sd combined-score sd" >clusters.scores
cat /dev/null >tt
foreach i (file.nam_clust[1-9]-scores file.nam_clust[1-9][0-9]-scores file.nam_clust[1-9][0-9][0-9]-scores)
  echo $i |sed -e 's/-scores//' -e 's/file\.nam_//' | awk '{printf "%-10s ",$1}' >>tt
  wc -l `echo $i| sed s/-scores//` |awk '{printf "%4d ",$1}'>>tt
  grep pdb $i | awk '{print $2}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %5.1f ",$1,$2}'>>tt
  grep pdb $i | awk '{print $3}' |$HADDOCKTOOLS/average.perl | awk '{printf "%7.1f %5.1f ",$1,$2}'>>tt
  grep pdb $i | awk '{print $4}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.0f %3.0f ",$1,$2}'>>tt
  grep pdb $i | awk '{print $5}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i | awk '{print $6}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i | awk '{print $7}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i | awk '{print $8}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i | awk '{print $9}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %5.1f ",$1,$2}'>>tt
  grep pdb $i | awk '{print $10}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.2f %5.2f ",$1,$2}'>>tt
  grep pdb $i | awk '{print $11}' |$HADDOCKTOOLS/average.perl | awk '{printf "%7.1f %7.1f\n",$1,$2}' >>tt
end
sort -g -k5 tt >>clusters.scores

if ($#argv > 0) then
echo "#cluster size rmsd_emin sd haddock-score sd bsa sd Edesolv sd fc-elec sd fc-desol sd Dfire sd probe sd probe/A**2 sd combined-score sd" >clusters_best{$1}ranked.scores
cat /dev/null >tt
foreach i (file.nam_clust[1-9]-scores file.nam_clust[1-9][0-9]-scores file.nam_clust[1-9][0-9][0-9]-scores)
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $1}' >`echo $i |sed -e 's/-scores/_best'$1'ranked/'`
  echo $i |sed -e 's/-scores//' -e 's/file\.nam_//' | awk '{printf "%-10s ",$1}' >>tt
  wc -l `echo $i| sed s/-scores//` |awk '{printf "%4d ",$1}'>>tt
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $2}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %5.1f ",$1,$2}'>>tt
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $3}' |$HADDOCKTOOLS/average.perl | awk '{printf "%7.1f %5.1f ",$1,$2}'>>tt
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $4}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.0f %3.0f ",$1,$2}'>>tt
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $5}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $6}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $7}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $8}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $9}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %5.1f ",$1,$2}'>>tt
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $10}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.2f %5.2f ",$1,$2}'>>tt
  grep pdb $i |sort -g -k11 |head -$1 | awk '{print $11}' |$HADDOCKTOOLS/average.perl | awk '{printf "%7.1f %7.1f\n",$1,$2}' >>tt
end
sort -g -k5 tt >>'clusters_best'$1'ranked.scores'

echo "#cluster size rmsd_emin sd haddock-score sd bsa sd Edesolv sd fc-elec sd fc-desol sd Dfire sd probe sd probe/A**2 sd combined-score sd" >clusters_best{$1}.scores
cat /dev/null >tt
foreach i (file.nam_clust[1-9]-scores file.nam_clust[1-9][0-9]-scores file.nam_clust[1-9][0-9][0-9]-scores)
  echo $i |sed -e 's/-scores//' -e 's/file\.nam_//' | awk '{printf "%-10s ",$1}' >>tt
  wc -l `echo $i| sed -e 's/-scores//'` |awk '{printf "%4d ",$1}'>>tt
  grep pdb $i |head -$1 | awk '{print $2}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %5.1f ",$1,$2}'>>tt
  grep pdb $i |head -$1 | awk '{print $3}' |$HADDOCKTOOLS/average.perl | awk '{printf "%7.1f %5.1f ",$1,$2}'>>tt
  grep pdb $i |head -$1 | awk '{print $4}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.0f %3.0f ",$1,$2}'>>tt
  grep pdb $i |head -$1 | awk '{print $5}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i |head -$1 | awk '{print $6}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i |head -$1 | awk '{print $7}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i |head -$1 | awk '{print $8}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %4.1f ",$1,$2}'>>tt
  grep pdb $i |head -$1 | awk '{print $9}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.1f %5.1f ",$1,$2}'>>tt
  grep pdb $i |head -$1 | awk '{print $10}' |$HADDOCKTOOLS/average.perl | awk '{printf "%5.2f %5.2f ",$1,$2}'>>tt
  grep pdb $i |head -$1 | awk '{print $11}' |$HADDOCKTOOLS/average.perl | awk '{printf "%7.1f %7.1f\n",$1,$2}' >>tt
end
sort -g -k5 tt >>clusters_best{$1}.scores
\rm tt

endif

