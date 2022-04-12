#!/bin/bash

# execute this script from inside the `examples/` folder.

# HADDOCK main relative path
HM="../../haddock3main/examples"

diff_capri () {

diff $1 ${HM}/$1

}

# ANTIBODY - ANTIGEN
echo "docking-antibody-antigen/CDR-acc-cltsel-test-capri_ss"
diff_capri docking-antibody-antigen/run1-CDR-acc-cltsel-test/12_caprieval/capri_ss.tsv

echo "docking-antibody-antigen/CDR-acc-cltsel-test-capri_clt"
diff_capri docking-antibody-antigen/run1-CDR-acc-cltsel-test/12_caprieval/capri_clt.tsv

echo "docking-antibody-antigen/CDR-acc-test-capri_ss"
diff_capri docking-antibody-antigen/run1-CDR-acc-test/11_caprieval/capri_ss.tsv

echo "docking-antibody-antigen/CDR-acc-test-capri_clt"
diff_capri docking-antibody-antigen/run1-CDR-acc-test/11_caprieval/capri_clt.tsv

echo "docking-antibody-antigen/ranairCDR-cltsel-test-capri_ss"
diff_capri docking-antibody-antigen/run1-ranairCDR-cltsel-test/12_caprieval/capri_ss.tsv

echo "docking-antibody-antigen/ranairCDR-cltsel-test-capri_ss"
diff_capri docking-antibody-antigen/run1-ranairCDR-cltsel-test/12_caprieval/capri_clt.tsv

echo "docking-antibody-antigen/ranairCDR-test-capri_ss"
diff_capri docking-antibody-antigen/run1-ranairCDR-test/8_caprieval/capri_ss.tsv

# PROTEIN - DNA
echo "docking-protein-DNA/run1-mdref-test"
diff_capri docking-protein-DNA/run1-mdref-test/7_caprieval/capri_ss.tsv

echo "docking-protein-DNA/run1-test"
diff_capri docking-protein-DNA/run1-test/7_caprieval/capri_ss.tsv

# PROTEIN - HOMOTRIMER
echo "docking-protein-homotrimer/run1-test"
diff_capri docking-protein-homotrimer/run1-test/7_caprieval/capri_ss.tsv

# PROTEIN - LIGAND
echo "docking-protein-ligand-shape/run1-test"
diff_capri docking-protein-ligand-shape/run1-test/5_caprieval/capri_ss.tsv

echo "docking-protein-ligand/run1-test"
diff_capri docking-protein-ligand/run1-test/5_caprieval/capri_ss.tsv

# PROTEIN - PEPTIDE
echo "docking-protein-peptide/run1-mdref-test"
diff_capri docking-protein-peptide/run1-mdref-test/7_caprieval/capri_ss.tsv

echo "docking-protein-peptide/run1-test"
diff_capri docking-protein-peptide/run1-test/7_caprieval/capri_ss.tsv

# PROTEIN - PROTEIN
echo "docking-protein-protein/run1-test"
diff_capri docking-protein-protein/run1-test/7_caprieval/capri_ss.tsv

echo "docking-protein-protein/run1-mdref-test"
diff_capri docking-protein-protein/run1-mdref-test/7_caprieval/capri_ss.tsv

echo "docking-protein-protein/run1-clt-test"
diff_capri docking-protein-protein/run1-clt-test/9_caprieval/capri_ss.tsv

# Refine complex
echo "refine-complex/run1-clt-test"
diff_capri refine-complex/run1-test/2_caprieval/capri_ss.tsv
