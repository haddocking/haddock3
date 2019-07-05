#!/bin/csh
#
if ($#argv == 0) then
  echo "Note that you can specify a number to perform "
  echo "the analysis for the Nbest structures of a cluster"
  echo "... analysis of fraction of native contacts of clusters running..."
else
  echo "... analysis of fraction of native contacts for "$1" best structures of each cluster..."
endif

foreach i (file.nam_clust[1-9] file.nam_clust[1-9][0-9] file.nam_clust[1-9][0-9][0-9] file.nam_clust[1-9]_best{$1} file.nam_clust[1-9][0-9]_best{$1} file.nam_clust[1-9][0-9][0-9]_best{$1} )
  cat /dev/null >$i'-fnat'
  foreach j (`cat $i`)
    grep $j file.nam_fnat >>$i'-fnat'
  end
end

echo "#cluster <fnat> sd " >clusters.fnat
cat /dev/null >tt
foreach i (file.nam_clust[1-9]-fnat file.nam_clust[1-9][0-9]-fnat file.nam_clust[1-9][0-9][0-9]-fnat)
  echo $i |sed -e 's/-fnat//' -e 's/file\.nam_//' | awk '{printf "%-10s ",$1}' >>tt
  grep pdb $i | awk '{print $2}' |$HADDOCKTOOLS/average.perl | awk '{printf "%6.2f %6.2f\n",$1,$2}' >>tt
end
sort -r -g -k2 tt >>clusters.fnat
\rm tt

if ($#argv > 0) then
echo "#cluster size <fnat> sd " >clusters_best{$1}.fnat
cat /dev/null >tt
foreach i (file.nam_clust[1-9]-fnat file.nam_clust[1-9][0-9]-fnat file.nam_clust[1-9][0-9][0-9]-fnat)
  echo $i |sed -e 's/-fnat//' -e 's/file\.nam_//' | awk '{printf "%-10s ",$1}' >>tt
  grep pdb $i | head -$1 | awk '{print $2}' |$HADDOCKTOOLS/average.perl | awk '{printf "%6.2f %6.2f\n",$1,$2}' >>tt
end
sort -r -g -k2 tt >>clusters_best{$1}.fnat
\rm tt

endif

