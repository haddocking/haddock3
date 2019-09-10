#
import json
import os
import glob

for folder in glob.glob('*'):
	pdb_a = f'{folder}/run1-ti-10r-1r/data/sequence/{folder}_r_u.pdb'
	pdb_b = f'{folder}/run1-ti-10r-1r/data/sequence/{folder}_l_u.pdb'
	ambig = f'{folder}/run1-ti-10r-1r/data/distances/ambig.tbl'
	run = {
		"input": {
			"molecules": {
				"mol1": f"{folder}_r_u.pdb",
				"mol2": f"{folder}_l_u.pdb"
			},
			"restraints": {
				"ambig": "ambig.tbl"
			},
			"run": 1,
			"nproc": 24,
			"run_scheme": "parallel"
		},
		"topology_recipe": "generate-topology.cns",
		"it0": {
			"recipe": "it0.cns",
			"sampling": 1000
		}
	}

	ref_a = f'/home/abonvin/BM5-full/{folder}/{folder}_r_b-matched.pdb'
	ref_b = f'/home/abonvin/BM5-full/{folder}/{folder}_l_b-matched.pdb'

	os.system(f'cat {ref_a} {ref_b} > {folder}/{folder}_complex.pdb')

	os.system(f'cp {pdb_a} {folder}/')
	os.system(f'cp {pdb_b} {folder}/')

	os.system(f'grep -v HETATM {folder}/{folder}_r_u.pdb > {folder}/oo')
	os.system(f'mv {folder}/oo {folder}/{folder}_r_u.pdb')

	os.system(f'grep -v HETATM {folder}/{folder}_l_u.pdb > {folder}/oo')
	os.system(f'mv {folder}/oo {folder}/{folder}_l_u.pdb')

	os.system(f'pdb_chain -A {folder}/{folder}_r_u.pdb > {folder}/oo')
	os.system(f'mv {folder}/oo {folder}/{folder}_r_u.pdb')
	os.system(f'python ~/pdb-tools/pdbtools/pdb_chainxseg.py {folder}/{folder}_r_u.pdb > {folder}/oo')
	os.system(f'mv {folder}/oo {folder}/{folder}_r_u.pdb')

	os.system(f'pdb_chain -B {folder}/{folder}_l_u.pdb > {folder}/oo')
	os.system(f'mv {folder}/oo {folder}/{folder}_l_u.pdb')
	os.system(f'python ~/pdb-tools/pdbtools/pdb_chainxseg.py {folder}/{folder}_l_u.pdb > {folder}/oo')
	os.system(f'mv {folder}/oo {folder}/{folder}_l_u.pdb')

	os.system(f'cp {ambig} {folder}/')
	with open(f'{folder}/run.json', 'w') as f:
		json.dump(run, f)

	with open(f'{folder}/run.job', 'w') as f:
		tbw = '#PBS -S /bin/tcsh\n'
		tbw += '#PBS -q short\n'
		tbw += '#PBS -l nodes=1:ppn=48\n'
		tbw += '#PBS -j oe\n'
		tbw += f'#PBS -N H3-it0-{folder}\n'
		tbw += '#PBS -m n\n'
		tbw += '\n'
		tbw += 'setenv PYTHONPATH $PYTHONPATH\:/home/rodrigo/haddock3\n'
		tbw += f'cd /home/rodrigo/work/haddock3-benchmark/{folder}\n'
		tbw += '/home/rodrigo/miniconda3/bin/python /home/rodrigo/haddock3/haddock/run_haddock.py run.json >& haddock.out\n'
		tbw += '\n'
		tbw += 'foreach pdb ( run1/it0/*pdb )\n'
		tbw += '  /home/software/haddock/haddock2.4/tools/pdb_segid-to-chain $pdb > oo\n'
		tbw += '  mv oo $pdb\n'
		tbw += 'end\n'
		tbw += '\n'
		tbw += 'foreach out ( run1/it0/*out )\n'
		tbw += '  gzip $out\n'
		tbw += 'end\n'
		tbw += '\n'
		tbw += f'/home/rodrigo/miniconda3/bin/python /home/rodrigo/haddock3/tools/capri-evaluate.py {folder}_complex' \
			f'.pdb run1/it0/file.list it0.haddock3.dat\n'
		tbw += f'/home/rodrigo/miniconda3/bin/python /home/rodrigo/haddock3/tools/capri-evaluate.py {folder}_complex' \
			f'.pdb run1-ti-10r-1r/structures/it0/file.list it0.haddock24.dat\n'
		tbw += '\n'
		tbw += 'touch DONE\n'
		tbw += '\n'
		f.write(tbw)
	f.close()
	os.system(f'chmod +x {folder}/run.job')
