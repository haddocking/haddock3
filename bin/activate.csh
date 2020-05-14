#!/bin/csh
setenv HADDOCK3 /home/abonvin/haddock3
setenv CNS_EXE /opt/nmr/cns_solve_1.31-UU/mac-intel-darwin/bin/cns

if (${?PYTHONPATH}) then
  setenv PYTHONPATH ${PYTHONPATH}:${HADDOCK3}
else
  setenv PYTHONPATH ${HADDOCK3}
endif
setenv PATH ${PATH}:${HADDOCK3}/bin