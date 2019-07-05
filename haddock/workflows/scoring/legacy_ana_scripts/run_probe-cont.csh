#!/bin/csh
#
setenv PROBEDIR /home/software/science/molprobity/bin

#echo "#structure probe-score probe-score/A**2" >probe-score.lis
foreach i ($argv)
  echo $i |awk '{printf "%s ",$1}' >>probe-score.lis
  set found=0
  if (-e $i:r.probe) then
    set found=`wc -l $i:r.probe | awk '{print $1}'`
  endif
  echo $found
  if ($found == 0) then
  pdb_xsegchain $i | grep -v DAN |grep -v ANI | sed -e 's/TYP/TYR/' | sed -e '/\ HA/d' \
    -e '/\ HB/d' \
    -e '/\ HG/d' \
    -e '/\ HD/d' \
    -e '/\ HE[0-9]/d' \
    -e '/\ HE\ /d' \
    -e '/\ HH/d' \
    -e '/\ HZ/d' \
    -e '/\ HT/d' \
    -e '/\ HN/d' \
    -e '/\ H8/d' \
    -e '/\ H7/d' \
    -e '/\ H6/d' \
    -e '/\ H5/d' \
    -e '/\ H4/d' \
    -e '/\ H3/d' \
    -e '/\ H2/d' \
    -e '/\ H1/d' \
    -e '/\ 1H/d' \
    -e '/\ 2H/d' \
    -e '/\ 3H/d' \
    -e '/\ H\ /d'\
    -e 's/OT1/O\ \ /' \
    -e 's/OT2/OXT/' >$i:r.pdb_seg
  $PROBEDIR/reduce -build -Xplor -quiet $i:r.pdb_seg >$i:rH.pdb
  $PROBEDIR/probe -quiet -MC -count -both "chainA" "chainB" $i:rH.pdb >$i:r.probe
  $PROBEDIR/probe -quiet -MC -count -both "chainA" "chainC" $i:rH.pdb >>$i:r.probe
  $PROBEDIR/probe -quiet -MC -count -both "chainB" "chainC" $i:rH.pdb >>$i:r.probe
  endif
  grep grand $i:r.probe | grep -v pdb | awk '{print $5}' | $HADDOCKTOOLS/sum.perl | awk '{printf "%6.1f ",$1}' >>probe-score.lis
  grep grand $i:r.probe | grep -v pdb | awk '{print $6}' | $HADDOCKTOOLS/sum.perl | awk '{printf "%6.2f\n",$1}' >>probe-score.lis
  if ($found == 0) then
    \rm $i:r.pdb_seg $i:rH.pdb
  endif
end

