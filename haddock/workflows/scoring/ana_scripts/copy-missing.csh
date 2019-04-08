#!/bin/csh -f
#
set num=$1

set fileroot=`grep fileroot run.cns | head -1 | awk 'BEGIN {FS=("=")}{print $5}' | sed s/\"//g |sed s/\;//`
set it0file='structures/it0/'`head -$num structures/it0/file.nam |tail -1`
set it0wfile=`echo $it0file |sed s/\.pdb/_water\.pdbw/`
set it1file='structures/it1/'$fileroot'_'$num'.pdb'
set it1wfile=`echo $it1file |sed s/\.pdb/_water\.pdbw/`
if (! -e $it1file) then
  cp $it0file $it1file
  if ( -e $it0wfile) then 
    cp $it0wfile $it1wfile
    if (`wc -l $it0wfile|awk '{print $1}'` == 0) then
      echo END >>$it1wfile
    endif
  endif
endif

