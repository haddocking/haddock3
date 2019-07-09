#!/bin/csh
#
setenv DFIRE /home/software/science/dcomplex

\cp $DFIRE/dfire.lib .
\cp $DFIRE/hash.dat .

echo "#structure Dfire-Ebinding[kcal/mol] dfire-score" >dfire-score.lis
foreach i ($argv)
   cat /dev/null >dfire
   $HADDOCKTOOLS/pdb_segid-to-chain $i |grep ATOM |sed -e 's/NEP/HIS/'|sed -e 's/CSP/CYS/' |grep -v CA2 | grep -v DAN |grep -v 'O[123]P' >$i:r.pdb_seg
   $DFIRE/dcomplex <<_Eod_ |sed -e 's/kcal\/mol//' |awk '{print $2,$3}' >>dfire
     $i:r.pdb_seg
     1
     A
     1
     B
_Eod_
   $DFIRE/dcomplex <<_Eod_ |sed -e 's/kcal\/mol//' |awk '{print $2,$3}' >>dfire
     $i:r.pdb_seg
     1
     A
     1
     C
_Eod_
   $DFIRE/dcomplex <<_Eod_ |sed -e 's/kcal\/mol//' |awk '{print $2,$3}' >>dfire
     $i:r.pdb_seg
     1
     B
     1
     C
_Eod_
   echo $i | awk '{printf "%s ",$1}' >>dfire-score.lis
   awk '{print $1}' dfire |$HADDOCKTOOLS/sum.perl |awk '{printf "%10.3f ",$1}' >>dfire-score.lis
   awk '{print $2}' dfire |$HADDOCKTOOLS/sum.perl |awk '{printf "%10.3f\n",$1}' >>dfire-score.lis
   \rm $i:r.pdb_seg dfire >&/dev/null
end
\rm dfire.lib
\rm hash.dat

