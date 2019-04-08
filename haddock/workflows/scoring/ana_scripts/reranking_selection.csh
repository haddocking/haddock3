#!/bin/csh

if (-e RERANKING_SELECTION_DONE) then
  exit
endif

\mv file.nam file.nam_original
\mv file.list file.list_original
\mv file.cns file.cns_original
if (`grep pdb structures.ranking |wc -l |awk '{print $1}'` != 0) then
  grep pdb structures.ranking-unsorted >file.ranking-unsorted
  grep pdb structures.ranking >file.ranking
  awk '{print $1}' file.ranking > file.nam
  paste file.list_original file.ranking-unsorted | awk 'NF==15 {print $1,$2,$15,$4}' | sort -nk3> file.list
  echo OK > RERANKING_SELECTION_DONE
  echo "RERANKING selection procedure successful..." > /dev/stderr
else
  echo "RERANKING selection procedure failed..." > /dev/stderr
  \mv file.nam_original file.nam
  \mv file.list_original file.list
  \mv file.cns_original file.cns
endif
