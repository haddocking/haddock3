#!/bin/csh
#
setenv WDIR /data/capri/Capri47/Target160-scoring/ana_scripts

foreach i (file.nam)
#  $WDIR/i-rmsd_to_xray.csh `head -2500 $i`
#  $WDIR/l-rmsd_to_xray.csh `head -2500 $i`
  $WDIR/run_fastcontact.csh `head -2500 $i`
  $WDIR/run_dfire.csh `head -2500 $i`
  $WDIR/run_probe.csh `head -2500 $i`
  set i1=`wc -l file.nam`
  @ i1-=2500
  if ($i1 > 0) then
#    $WDIR/i-rmsd_to_xray-cont.csh `tail -$i1 $i`
#    $WDIR/l-rmsd_to_xray-cont.csh `tail -$i1 $i`
    $WDIR/run_fastcontact-cont.csh `tail -$i1 $i`
    $WDIR/run_dfire-cont.csh `tail -$i1 $i`
    $WDIR/run_probe-cont.csh `tail -$i1 $i`
  endif
  $WDIR/ana_structures.csh
  paste structures_haddock-sorted.stat fastcontact-score.lis dfire-score.lis probe-score.lis| head -1 |\
    awk '{print $1,$3,$2,$5,$6,$7,$8,$9,$21,$22,$23,$25,$26,$28,$31,$32}' >structures.scores
  paste structures_haddock-sorted.stat fastcontact-score.lis dfire-score.lis probe-score.lis| grep -v "#" |\
    awk '{print $1,$3,$2,$5,$6,$7,$8,$9,$21,$22,$23,$25,$26,$28,$31,$32}' >>structures.scores
  head -1 structures.scores | awk '{print $1,$2,$3,$9,$11,$12,$13,$14,$15,$16,"combined-score"}' > structures.ranking
  grep pdb structures.scores | awk '{if ($15<=0){sc=0.01}else{sc=$15};print $1,$2,$3,$9,$11,$12,$13,$14,$15,$16,($3+$12+(5*$13)+0*$14-0*$15)}'| sort -g -k11>>structures.ranking
  head -1 structures.scores | awk '{print $1,$2,$3,$9,$11,$12,$13,$14,$15,$16,"combined-score"}' > structures.ranking-unsorted
  grep pdb structures.scores | awk '{if ($15<=0){sc=0.01}else{sc=$15};print $1,$2,$3,$9,$11,$12,$13,$14,$15,$16,($3+$12+(5*$13)+0*$14-0*$15)}' >> structures.ranking-unsorted
end
