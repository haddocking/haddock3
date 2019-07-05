#!/bin/csh
#
# Check for awk or gawk
#
set found=`which gawk |grep -v found |wc -l`
if ($found == 0) then
  set found=`which awk |grep -v found |wc -l`
  if ($found == 0) then
    echo 'awk or gawk not found'
    echo '==> no analysis possible'
    goto exit
  else
    set AWK=awk
  endif
else
  set AWK=gawk
endif
#
# Define the location of profit
#
if ( `printenv |grep PROFIT | wc -l` == 0) then
  set found=`which profit |wc -l`
  if ($found == 0) then
     echo 'PROFIT environment variable not defined'
     echo '==> no rmsd calculations '
  else
     setenv PROFIT `which profit`
  endif
endif
#
# extract HADDOCK score
#
foreach iclu (file.nam)
  echo "#struc haddock-score" >$iclu'_haddock-score'
  foreach i (`head -5000 $iclu`)
    echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_haddock-score'
    grep $i file.list | $AWK '{print $3}'>>$iclu'_haddock-score'
  end
  set i1=`wc -l file.nam`
  @ i1-=5000
  if ($i1 > 0) then
    foreach i (`tail -$i1 $iclu`)
      echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_haddock-score'
      grep $i file.list | $AWK '{print $3}'>>$iclu'_haddock-score'
    end
  endif
end
#
# Extract buried surface area
#
foreach iclu (file.nam)
  echo "#struc bsa" >$iclu'_bsa'
  foreach i (`head -5000 $iclu`)
    echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_bsa'
    head -30 $i | grep -i buried | $AWK '{print $5}'>>$iclu'_bsa'
  end
  set i1=`wc -l file.nam`
  @ i1-=5000
  if ($i1 > 0) then
    foreach i (`tail -$i1 $iclu`)
      echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_bsa'
      head -30 $i | grep -i buried | $AWK '{print $5}'>>$iclu'_bsa'
    end
  endif
end
#
# Extract energies
#
foreach iclu (file.nam)
  echo "#struc Einter  Enb Evdw+0.1Eelec Evdw  Eelec Eair Ecdih Ecoup Esani Evean Edani" >$iclu'_ener'
  foreach i (`head -5000 $iclu`)
    echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_ener'
    head -30 $i | grep -i energies | sed -e's/\,//g' | $AWK '{if(NF==15){print $3,$8+$9,$8+0.1*$9,$8,$9,$10,$11,$12,$13,$14,$15} else {print $3,$8+$9,$8+0.1*$9,$8,$9,$10,$11,$12,$13,$14,0.0}}'>>$iclu'_ener'
  end
  set i1=`wc -l file.nam`
  @ i1-=5000
  if ($i1 > 0) then
    foreach i (`tail -$i1 $iclu`)
      echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_ener'
      head -30 $i | grep -i energies | sed -e's/\,//g' | $AWK '{if(NF==15){print $3,$8+$9,$8+0.1*$9,$8,$9,$10,$11,$12,$13,$14,$15} else {print $3,$8+$9,$8+0.1*$9,$8,$9,$10,$11,$12,$13,$14,0.0}}'>>$iclu'_ener'
    end
  endif
end
#
# extract violations
#
foreach iclu (file.nam)
  echo "#struc #NOEviol #Dihedviol #Jviol #Saniviol #veanviol #Daniviol" >$iclu'_viol'
  foreach i (`head -5000 $iclu`)
    echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_viol'
    head -30 $i | grep -i " violation" | sed -e's/\,//g' | $AWK '{if (NF==8){print $3,$4,$5,$6,$7,$8} else {print $3,$4,$5,$6,$7,0}}'>>$iclu'_viol'
  end
  set i1=`wc -l file.nam`
  @ i1-=5000
  if ($i1 > 0) then
    foreach i (`tail -$i1 $iclu`)
      echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_viol'
      head -30 $i | grep -i " violation" | sed -e's/\,//g' | $AWK '{if (NF==8){print $3,$4,$5,$6,$7,$8} else {print $3,$4,$5,$6,$7,0}}'>>$iclu'_viol'
    end
  endif
end
#
# extract binding energy
#
foreach iclu (file.nam)
  echo "#struc dH" >$iclu'_dH'
  foreach i (`head -5000 $iclu`)
    echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_dH'
    head -30 $i | grep -i " Binding" | $AWK '{print $4}'>>$iclu'_dH'
  end
  set i1=`wc -l file.nam`
  @ i1-=5000
  if ($i1 > 0) then
    foreach i (`tail -$i1 $iclu`)
      echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_dH'
      head -30 $i | grep -i " Binding" | $AWK '{print $4}'>>$iclu'_dH'
    end
  endif
end
#
# extract desolvation energy
#
foreach iclu (file.nam)
  echo "#struc Edesolv" >$iclu'_Edesolv'
  foreach i (`head -5000 $iclu`)
    echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_Edesolv'
    head -30 $i | grep -i " desolv" | $AWK '{print $4}'>>$iclu'_Edesolv'
  end
  set i1=`wc -l file.nam`
  @ i1-=5000
  if ($i1 > 0) then
    foreach i (`tail -$i1 $iclu`)
      echo $i |$AWK '{printf "%s ",$1}'>>$iclu'_Edesolv'
      head -30 $i | grep -i " desolv" | $AWK '{print $4}'>>$iclu'_Edesolv'
    end
  endif
end
#
# RMSD calculations only if PROFITDIR is defined
#
if (-e $PROFIT) then

if (-e i-RMSD.dat) then
  \cp i-RMSD.dat file.nam_rmsd
else
#
# calculate rmsd from lowest energy model (bakcbone -start and end residues)
#
set refe=`head -1 file.nam`
set atoms='CA'

cat /dev/null > rmsd_best.disp
foreach iclu (file.nam)
  foreach i (`head -5000 $iclu`)
    echo $i >>rmsd_best.disp
    $PROFIT <<_Eod_ |grep RMS |tail -1>>rmsd_best.disp
      refe $refe
      mobi $i
      ignore missing
      atom $atoms
      zone A*
      fit
      rzone B*
      quit
_Eod_
  end
  set i1=`wc -l file.nam`
  @ i1-=5000
  if ($i1 > 0) then
    foreach i (`tail -$i1 $iclu`)
      echo $i >>rmsd_best.disp
      $PROFIT <<_Eod_ |grep RMS >>rmsd_best.disp
        refe $refe
        mobi $i
	ignore missing
	zone A*
        atom $atoms
        fit
	rzone B*
        quit
_Eod_
    end
  endif
  echo ' ' >>rmsd_best.disp
  echo "#structure rmsd_all" >$iclu'_rmsd'
  $AWK '{if ($1 == "RMS:") {printf "%8.3f ",$2} else {printf "\n %s ",$1}}' rmsd_best.disp |grep pdb >>$iclu'_rmsd'
  echo ' ' >>$iclu'_rmsd'
  awk '{if (NF == 1) {$2=999}; print $1,$2}' $iclu'_rmsd' >ttt
  \mv ttt $iclu'_rmsd'
  \rm rmsd_best.disp
end
endif
#
# Combine results in one file
# 
foreach iclu (file.nam)
  paste $iclu'_haddock-score' $iclu'_rmsd' $iclu'_ener' $iclu'_viol' $iclu'_bsa' $iclu'_dH' $iclu'_Edesolv'\
 | $AWK '{print $1,$2,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$18,$19,$20,$21,$22,$23,$25,$27,$29}' >structures_haddock-sorted.stat
end
grep "#struc" structures_haddock-sorted.stat > structures_rmsd-sorted.stat
sort -n -k3 structures_haddock-sorted.stat |grep pdb >> structures_rmsd-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_ene-sorted.stat
sort -n -k4 structures_haddock-sorted.stat |grep pdb >> structures_ene-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_nb-sorted.stat
sort -n -k5 structures_haddock-sorted.stat |grep pdb >> structures_nb-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_nbw-sorted.stat
sort -n -k6 structures_haddock-sorted.stat |grep pdb >> structures_nbw-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_air-sorted.stat
sort -n -k9 structures_haddock-sorted.stat |grep pdb >> structures_air-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_airviol-sorted.stat
sort -n -k15 structures_haddock-sorted.stat |grep pdb >> structures_airviol-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_bsa-sorted.stat
sort -r -n -k21 structures_haddock-sorted.stat |grep pdb >> structures_bsa-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_dH-sorted.stat
sort -n -k22 structures_haddock-sorted.stat |grep pdb >> structures_dH-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_Edesolv-sorted.stat
sort -n -k23 structures_haddock-sorted.stat |grep pdb >> structures_Edesolv-sorted.stat
#
# else no RMSD calculation
#
else
#
# Combine results in one file
# 
foreach iclu (file.nam)
  paste $iclu'_haddock-score' $iclu'_ener' $iclu'_viol' $iclu'_bsa'  $iclu'_dH' $iclu'_Edesolv'\
 | $AWK '{print $1,$2,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16,$17,$18,$19,$20,$21,$23,$25,$27}' >structures_haddock-sorted.stat
end
grep "#struc" structures_haddock-sorted.stat > structures_ene-sorted.stat
sort -n -k3 structures_haddock-sorted.stat |grep pdb >> structures_ene-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_nb-sorted.stat
sort -n -k4 structures_haddock-sorted.stat |grep pdb >> structures_nb-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_nbw-sorted.stat
sort -n -k5 structures_haddock-sorted.stat |grep pdb >> structures_nbw-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_air-sorted.stat
sort -n -k8 structures_haddock-sorted.stat |grep pdb >> structures_air-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_airviol-sorted.stat
sort -n -k14 structures_haddock-sorted.stat |grep pdb >> structures_airviol-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_bsa-sorted.stat
sort -r -n -k20 structures_haddock-sorted.stat |grep pdb >> structures_bsa-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_dH-sorted.stat
sort -n -k21 structures_haddock-sorted.stat |grep pdb >> structures_dH-sorted.stat
grep "#struc" structures_haddock-sorted.stat > structures_Edesolv-sorted.stat
sort -n -k22 structures_haddock-sorted.stat |grep pdb >> structures_Edesolv-sorted.stat

endif
