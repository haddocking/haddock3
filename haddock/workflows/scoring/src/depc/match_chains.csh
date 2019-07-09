#!/bin/csh
#syntax: match_chains.csh <list of pdb files>
#make sure that the pdb files are in a remote directory

setenv HADDOCKSCRIPTS /home/software/haddock/haddock_scripts/PDB-matching
source ~abonvin/haddock2.4/haddock_configure.csh

mkdir pdb-matched
foreach i(`cat $1`)
  $HADDOCKSCRIPTS/pdb_chain-segid $i > pdb-matched/$i:t
end

set refe = `head -1 $1`
set refe = $refe:t

cd pdb-matched
#if (0) then

rm -f chain_convert_?.list
$HADDOCKSCRIPTS/match_chains.py $refe $refe | awk '{system("cat /dev/null > chain_convert_" $1 ".list")}'
foreach i(`cat ../$1`)
 set ii = $i:t
 $HADDOCKSCRIPTS/match_chains.py $refe $ii | awk -v file=$ii '$2 == "None" {print "Error: " file " has no match for chain $1 "; next} {tar = "chain_convert_" $1 ".list"; print file, $2 >> tar}'
end

#endif

#if (0) then
foreach c(chain_convert_?.list)
  set chain = `echo $c:r | sed 's/chain_convert_//'`
  \cp $refe common_$chain.pdb  
  foreach i(`sed 's/ /+/g' $c`)
    set ii = `echo $i | sed 's/+/ /g' | awk '{print $1}'`
    set ichain = `echo $i | sed 's/+/ /g' | awk '{print $2}'`
    $HADDOCKSCRIPTS/pdb-pdbalign common_$chain.pdb $chain $ii $ichain |& awk '$1 != "Warning:"' | awk '$1 == "ATOM" && substr($0,22,1) != "X"' > tmp
#    \cp common_$chain.pdb common_$chain.pdb_$ii
    \mv tmp common_$chain.pdb
  end
end
#endif


foreach c(chain_convert_?.list)
  set chain = `echo $c:r | sed 's/chain_convert_//'`
  foreach i(`sed 's/ /+/g' $c`)
    set ii = `echo $i | sed 's/+/ /g' | awk '{print $1}'`
    set ichain = `echo $i | sed 's/+/ /g' | awk '{print $2}'`  
    $HADDOCKSCRIPTS/pdb-pdbalign common_$chain.pdb $chain $ii $ichain |& awk '$1 != "Warning:"' | awk '$1 == "ATOM" && substr($0,22,1) != "X"' > {$ii}_$chain
  end
end

foreach i(`cat ../$1`)
  cat {$i:t}_? > ttt
  $HADDOCKTOOLS/pdb_chain-to-segid ttt > $i:t
  echo END >>$i:t
  rm -f {$i:t}_? ttt
end

cd ../
