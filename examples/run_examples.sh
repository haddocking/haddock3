#!/bin/bash

# Script used to run all examples at once.
# Does not parse examples automatically. The exact examples that are run
# are defined here in the script. Run this script from inside the `examples`
# folder.
#
# USAGE:
#
# ./run_examples.sh 0
# ./run_examples.sh 1


argmsg=$(cat <<-END
Please give 0 to run all examples regardless of errors.
Or, give 1, to run all examples but stop if an error happens.
END
)


if [ $# -eq 0 ]; then
    echo "No arguments provided."
    echo ${argmsg}
    exit 1
fi

if [[ ! "0 1" =~ $1 ]]; then
    echo  "Input argument not recognized."
    echo ${argmsg}
    exit 1
fi


e1 () {
cd docking-protein-DNA
rm -r run1
echo "************** RUNNING PROTEIN-PROTEIN-DNA *************"
haddock3 docking-protein-DNA.cfg
}

e2 () {
cd docking-protein-DNA
rm -r run1-mdref
echo "************** RUNNING PROTEIN-PROTEIN-DNA-MDREF *************"
haddock3 docking-protein-DNA-mdref.cfg
}

e3 () {
cd docking-protein-homotrimer
rm -r run1
echo "************** RUNNING PROTEIN-HOMOTRIMER *************"
haddock3 docking-protein-homotrimer.cfg
}

e4 () {
cd docking-protein-ligand-shape
rm -r run1
echo "************** RUNNING PROTEIN-LIGAND-SHAPE *************"
haddock3 docking-protein-ligand-shape.cfg
}

e5 () {
cd docking-protein-ligand
rm -r run1
echo "************** RUNNING PROTEIN-LIGAND *************"
haddock3 docking-protein-ligand.cfg
}

e6 () {
cd docking-protein-peptide
rm -r run1
echo "************** RUNNING PROTEIN-PEPTIDE *************"
haddock3 docking-protein-peptide.cfg
}

e7 () {
cd docking-protein-protein
rm -r run1
echo "************** RUNNING PROTEIN-PROTEIN *************"
haddock3 docking-protein-protein.cfg
}

e8 () {
cd docking-protein-protein
echo "************** RUNNING PROTEIN-PROTEIN-MDREF *************"
rm -r run1-mdref
haddock3 docking-protein-protein-mdref.cfg
}

e9 () {
cd refine-complex
rm -r run1
echo "************** RUNNING REFINE-COMPLEX *************"
haddock3 refine-complex.cfg
}

e10 () {
cd scoring
rm -r run1
echo "************** RUNNING SCORING *************"
haddock3 scoring.cfg
}

examples=( e1 e2 e3 e4 e5 e6 e7 e8 e9 e10 )

for example in ${examples[@]}; do
    echo ""
    ${example}
    exitcodefrompreviouscommand=$?  # ;-)
    if [ ${exitcodefrompreviouscommand} -ne 0 ]; then
        if [ ${1} -eq 1 ]; then
            break
        elif [ ${1} -eq 0 ]; then
            cd ..
            continue
        fi
    fi
    cd ..
done
exit ${exitcodefrompreviouscommand}
