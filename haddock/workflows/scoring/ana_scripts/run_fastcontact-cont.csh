#!/bin/csh
#
setenv FASTCONTACT /home/software/science/fastcontact

foreach i ($argv)
   $HADDOCKTOOLS/pdb_segid-to-chain $i >tmp
   pdb_selxseg -A tmp |grep -v CA2 |sed -e 's/\ HN\ /\ H\ \ /' -e 's/CSP/CYS/'|grep ATOM >a.pdb
   pdb_selxseg -B tmp |grep -v CA2 |sed -e 's/\ HN\ /\ H\ \ /' -e 's/CSP/CYS/'|grep ATOM >b.pdb
   pdb_selxseg -C tmp |grep -v CA2 |sed -e 's/\ HN\ /\ H\ \ /' -e 's/CSP/CYS/'|grep ATOM >c.pdb
   echo $i |awk '{printf "%s ",$1}' >>fastcontact-score.lis
   $FASTCONTACT/fastcontact a.pdb b.pdb |head -2 |tail -1 >tmp
   cat tmp |awk '{print $1}' >felec
   cat tmp |awk '{print $2}' >fdesol
   $FASTCONTACT/fastcontact a.pdb c.pdb |head -2 |tail -1 >tmp
   cat tmp |awk '{print $1}' >>felec
   cat tmp |awk '{print $2}' >>fdesol
   $FASTCONTACT/fastcontact b.pdb c.pdb |head -2 |tail -1 >tmp
   cat tmp |awk '{print $1}' >>felec
   cat tmp |awk '{print $2}' >>fdesol
   $HADDOCKTOOLS/sum.perl felec | awk '{printf "%10.3f ",$1}' >>fastcontact-score.lis
   $HADDOCKTOOLS/sum.perl fdesol | awk '{printf "%10.3f\n",$1}' >>fastcontact-score.lis
   \rm a.pdb b.pdb c.pdb d.pdb tmp fort.* felec fdesol >&/dev/null
end

