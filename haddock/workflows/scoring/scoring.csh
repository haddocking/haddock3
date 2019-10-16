#!/usr/bin/env tcsh
#PBS -N exp-scoring
#PBS -q short
#PBS -l nodes=1:ppn=48
#PBS -S /bin/tcsh

echo $PYTHONPATH
echo `pwd`

setenv PYTHONPATH $PYTHONPATH\:/home/rodrigo/haddock3
setenv PYTHONPATH $PYTHONPATH\:/home/rodrigo/haddock3/haddock/workflows/scoring
cd /home/rodrigo/work/scoring
/home/rodrigo/miniconda3/bin/python /home/rodrigo/haddock3/haddock/workflows/scoring/run_scoring.py /home/rodrigo/work/scoring/scoring.json
