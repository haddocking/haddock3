#!/bin/csh
#
source ../../haddock_configure.csh
echo "=========================================================="
echo "=========================================================="
echo " RUNNING NOW E2A-HPR PROTEIN-PROTEIN DOCKING EXAMPLE"
echo "=========================================================="
echo "=========================================================="
haddock2.4 >&/dev/null
cd run1
patch -p0 -i ../run.cns.patch >&/dev/null
haddock2.4 >&haddock.out
cd ..
./ana_scripts/run_all.csh run1 >&/dev/null
../results-stats.csh run1
echo "=========================================================="
echo "=========================================================="
echo " E2A-HPR PROTEIN-PROTEIN DOCKING EXAMPLE COMPLETED"
echo "=========================================================="
echo "=========================================================="
