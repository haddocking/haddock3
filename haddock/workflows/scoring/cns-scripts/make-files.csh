#!/bin/csh
#
cat /dev/null >t1
cat /dev/null >t2
foreach i ($argv)
  echo $i |awk '{printf "%s ",$1}' >>t1
  grep energies $i | sed s/\,/\ /g | awk '{printf "%f %f\n",$8,$9}' >>t1
  grep Desol $i | sed s/\,/\ /g | awk '{print $4}' >>t2
  grep Esym $i | sed s/\,/\ /g | awk '{print $4}' >>t3
end
paste t1 t2 t3 | awk '{printf "%s%s %f%s\n", $1," { ",($2+0.2*$3+$4+0.1*$5)," }"}' |sort -n +2 >file.list
awk '{print $1}' file.list >file.nam
\rm t1 t2 t3
