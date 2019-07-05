#!/bin/csh
#
cat /dev/null >t1
cat /dev/null >t2
cat /dev/null >t3
cat /dev/null >t4
foreach i ($argv)
  echo $i |awk '{printf "%s ",$1}' >>t1
  grep energies $i | sed s/\,/\ /g | awk '{printf "%f %f\n",$8,$9}' >>t1
  grep buried $i | sed s/\,/\ /g | awk '{print $5}' >>t2
  grep Desol $i | sed s/\,/\ /g | awk '{print $4}' >>t3
  grep Esym $i | sed s/\,/\ /g | awk '{print $4}' >>t4
end
paste t1 t2 t3 t4 | awk '{printf "%s%s %f%s\n", $1," { ",(0.01*$2+$3-0.01*$4+$5+0.1*$6)," }"}' |sort -n +2 >file.list
awk '{print $1}' file.list >file.nam
\rm t1 t2 t3 t4
