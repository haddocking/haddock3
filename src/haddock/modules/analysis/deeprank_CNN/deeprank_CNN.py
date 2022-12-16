from mpi4py import MPI
import os
import json
import requests
import pandas as pd
import numpy as np
import torch
import torch.nn.functional as F
import h5py
from Bio import pairwise2
import deeprank
from deeprank.generate import *
from deeprank.learn import NeuralNet
from haddock.modules.analysis.deeprank_CNN.model_280619 import cnn_class

from requests.packages.urllib3.exceptions import InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)


#used to rank docking poses

class CNN_score():
    def __init__(self, pdb_source, chain1, chain2, output_dir):
        self.pdb_source = pdb_source
        self.chain1 = chain1
        self.chain2 = chain2
        self.output_dir = output_dir
        self.pssm_path = os.path.join(pdb_source, 'pssm')
        
        for chainID in [self.chain1, self.chain2]:
            self.write_pssm(chainID)
        for chainID in [self.chain1, self.chain2]:
            fpssm = [file for file in os.listdir(self.pssm_path) if f'{chainID}.'in file][0]
            pdbs = [file for file in os.listdir(self.pdb_source) if '.pdb' in file]
            for fpdb in pdbs:
                self.write_mapped_pssm_pdb(fpssm, fpdb, chainID, self.pssm_path)
    
    def write_pssm(self, chainID):
        '''
        query pssm database
        '''
        pdb_code = [f for f in os.listdir(self.pdb_source) if '.pdb' in f][0].split('_')[0]
        url = f'https://3dcons.cnb.csic.es/pssm_json/{pdb_code}/{chainID}'
        pssm_dict = requests.get(url, verify=False).json()
        if pssm_dict is not None:
            header = 'pdbresi pdbresn seqresi seqresn A R N D C Q E G H I L K M F P S T W Y V\n'
            for i in pssm_dict:
                pdbresi = str(i['res_id'])
                pdbresn = i['aa']
                seqresi = str(i['index'])
                seqresn = i['aa']
                pssm = ' '.join(i['iter']['3']['pssm'])
                line = pdbresi+ ' ' + pdbresn + ' '+ seqresi + ' ' + seqresn + ' '+ pssm + '\n'
                header += line
            if not os.path.isdir(self.pssm_path):
                os.mkdir(self.pssm_path)
            with open(os.path.join(self.pssm_path, f'{pdb_code}_{chainID}.pssm'), 'w') as f:
                f.write(header)
        else:
            raise ValueError(
            "Please enter a valid pdb code and chainID")

    def get_pssm(self, fpssm):
        '''
        read pssm file
        '''
        rule = tuple([str(i) for i in range(10)])
        pssm = []
        with open(os.path.join(self.pssm_path, fpssm), "r") as f:
            for line in f.readlines():
                line_raw = line
                line = line.strip()
                # only select lines that contain pssm values
                if line.startswith(rule):
                    if len(line.split()) == 24:
                        pssm.append(line.split())
                    else:
                        raise ValueError(
                            "Wrong format of the following line in PSSM file {}:\n{}".format(fpssm, line_raw))
        return pssm

    def seq_from_pssm(self, fpssm):
        '''
        get sequence from pssm
        '''
        seqs = []
        res_index = []
        with open(os.path.join(self.pssm_path, fpssm), "r") as f:
            for line in f.readlines()[1:]:
                seqs.append(line.split()[1])
                res_index.append(line.split()[0])
        return seqs, res_index

    def seq_from_pdb(self, fpdb, chainID):
        '''
        get sequence from pdb
        '''
        res_codes = [
            # 20 canonical amino acids
            ('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
            ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
            ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
            ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
            ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
            # Non-canonical amino acids
            ('ASX', 'B'), ('SEC', 'U'), ('GLX', 'Z'),
            # ('MSE', 'M'), ('SOC', 'C'),
            # Canonical xNA
            ('  U', 'U'), ('  A', 'A'), ('  G', 'G'), ('  C', 'C'),
            ('  T', 'T'),
        ]

        three_to_one = dict(res_codes)
        _records = set(['ATOM  ', 'HETATM'])

        chainID = chainID.upper()
        sequence = []
        resID = []
        chains = set()
        read = set()
        with open(os.path.join(self.pdb_source, fpdb), "r") as f:
            for line in f:
                line = line.strip()
                if line[0:6] in _records:
                    resn = line[17:20]
                    chain = line[21]
                    resi = line[22:26].replace(' ', '')
                    icode = line[26]
                    r_uid = (resn, chain, resi, icode)
                    chains.add(chain)
                    if chain == chainID:
                        if r_uid not in read:
                            read.add(r_uid)
                        else:
                            continue
                        aa_resn = three_to_one.get(resn, 'X')
                        sequence.append(aa_resn)
                        resID.append(resi)
            if chainID not in chains:
                raise ValueError(
                    "Chain `{}` NOT exist in PDB file '{}'".format(chainID, fpdb))
        return sequence, resID

    
    def get_aligned_sequences(self, pdb_seq, pssm_seq):
        '''
        align protein sequence
        '''
        ali = pairwise2.align.globalxs(pdb_seq, pssm_seq, -2, -1)
        seq1_ali = np.array([i for i in ali[0][0]])
        seq2_ali = np.array([i for i in ali[0][1]])

        return seq1_ali, seq2_ali

    def write_mapped_pssm_pdb(self, fpssm, fpdb, chainID, pssm_path):
        '''
        solve protein numbering problem
        adapted from pssmgen package
        '''
        #pssmname = os.path.basename(fpssm)
        pdbname = os.path.basename(fpdb)
        pssm_seq, _= self.seq_from_pssm(fpssm)
        pdb_seq, pdb_resn = self.seq_from_pdb(fpdb, chainID)
        pssm_seq= ''.join(pssm_seq)
        pdb_seq = ''.join(pdb_seq)
        pdb_seq_align, pssm_seq_align = self.get_aligned_sequences(pdb_seq, pssm_seq)
        index_match = pdb_seq_align == pssm_seq_align
        index_mismatch = np.logical_not(index_match)
        seqlen = len(pdb_seq_align)

        gap_seq = np.array(["-"] * seqlen)
        resX_seq = np.array(["X"] * seqlen)
        index_gappdb = gap_seq == pdb_seq_align
        index_resXpdb = resX_seq == pdb_seq_align
        index_gappssm = gap_seq == pssm_seq_align
        index_resXpssm = resX_seq == pssm_seq_align
        # get index of normal residues (not gap, not res X) for each sequence
        index_norm_pdb = np.logical_not(np.logical_or(index_gappdb, index_resXpdb))
        index_norm_pssm = np.logical_not(np.logical_or(index_gappssm, index_resXpssm))
        # get index of normal residues for both sequences
        index_norm = np.logical_and(index_norm_pdb, index_norm_pssm)
        # get index of mutated normal residues
        index_mut = np.logical_and(index_mismatch, index_norm)

        if len(set(index_mut)) > 1:
            mut_seq = []
            for i in index_mut:
                if i:
                    mut_seq.append("^")
                else:
                    mut_seq.append("_")

        index_norm_nogappssm = index_norm[np.logical_not(index_gappssm)]
        pssm = np.array(self.get_pssm(fpssm))
        pssm_norm = pssm[index_norm_nogappssm]

        # add the residue number and name of PDB file to the mapped pssm
        # for pssm content, only keep the scoring matrix and information content
        header = ["pdbresi", "pdbresn", "seqresi", "seqresn", "A", "R", "N",
                "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P",
                "S", "T", "W", "Y", "V", "IC"]
        header = np.transpose(np.array([[i] for i in header]))

        pdb_resn = [[i] for i in pdb_resn]
        pdb_seq = [[i] for i in pdb_seq]
        index_norm_nogappdb = index_norm[np.logical_not(index_gappdb)]
        resi_pdb = np.array(pdb_resn)[index_norm_nogappdb]
        resn_pdb = np.array(pdb_seq)[index_norm_nogappdb]
        pssm_out = np.concatenate((resi_pdb, resn_pdb, pssm_norm[:,:2], pssm_norm[:, 4:24], pssm_norm[:, -2:-1]), axis=1)
        pssm_out = np.concatenate((header, pssm_out))

        # write mapped pssm to file which is named with input PDB file name, chain ID and ".pdb.pssm"
        fopssm = os.path.join(self.pssm_path, os.path.splitext(pdbname)[0] + "." + chainID.upper() + ".pssm")

        with open(fopssm, "w") as f:
            for i in pssm_out:
                tmp1 = ["{:>7s}".format(j) for j in i[:4]]
                tmp2 = ["{:>4s}".format(j) for j in i[4:]]
                f.write(" ".join(tmp1+tmp2) + "\n")
        
    '''
    #old local pssm generation script
    #produce the exact same format as in deeprank
    def generate_pssm(self):
        gen = PSSM(work_dir=self.pdb_source)
        gen.configure(blast_exe='/trinity/login/xxu/software/ncbi-blast-2.13.0+/bin/psiblast',
            database='/trinity/login/xxu/data/DBs/nr',
            num_threads=4, evalue=0.0001, comp_based_stats='T',
            max_target_seqs=2000, num_iterations=3, outfmt=7,
            save_each_pssm=True, save_pssm_after_last_round=True)
        gen.get_fasta(pdb_dir='', chain=(self.chain1,self.chain2), out_dir='fasta')
        gen.get_pssm(fasta_dir='fasta', out_dir='pssm_raw', run=True)
        gen.map_pssm(pssm_dir='pssm_raw', pdb_dir='', out_dir='pssm', chain=(self.chain1,self.chain2))
        gen.get_mapped_pdb(pdbpssm_dir='pssm', pdb_dir='', pdbnonmatch_dir='pdb_nonmatch')
    '''

    def create_database(self):
        '''
        create deeprank database
        '''
        comm = MPI.COMM_WORLD
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
        self.hdf_dir = os.path.join(self.output_dir, 'output.hdf5')
        database = DataGenerator(pdb_source= self.pdb_source, #path to the models  
                         pssm_source=self.pssm_path, #path to the pssm data
                         data_augmentation = None,
                         chain1=self.chain1, chain2=self.chain2,
                         compute_features = ['deeprank.features.AtomicFeature', 'deeprank.features.FullPSSM','deeprank.features.PSSM_IC'
                         , 'deeprank.features.BSA', 'deeprank.features.ResidueDensity'],
                         hdf5=self.hdf_dir,mpi_comm=comm)
        database.create_database(prog_bar=True)

        grid_info = {
            'number_of_points' : [30,30,30],
            'resolution' : [1.,1.,1.],
            'atomic_densities': {'C': 1.7, 'N': 1.55, 'O': 1.52, 'S': 1.8}
            }
        database.map_features(grid_info,try_sparse=True, time=False, prog_bar=True)

    def test_CNN(self):
        '''
        apply pre-trained model
        '''
        model_data = 'best_train_model.pt'
        path = '/trinity/login/xxu/software/haddock3/src/haddock/modules/analysis/deeprank_CNN'
        model_data = os.path.join(path, model_data)
        database = self.hdf_dir
        model = NeuralNet(database, 
                         cnn_class, 
                         task='class', 
                         pretrained_model=model_data, 
                         save_hitrate=False,
                         outdir=self.output_dir)
        model.test(hdf5='prediction.hdf5')

    def analysis_result(self):
        path = os.path.join(self.output_dir, 'prediction.hdf5')
        f = h5py.File(path)
        out_class = f['epoch_0000']['test']['outputs'][()]
        out = F.softmax(torch.FloatTensor(out_class),dim=1).data.numpy()[:, 1]
        ind_sort = np.argsort(out)[::-1]
        out = out[ind_sort]
        probility = F.softmax(torch.FloatTensor(out_class), dim=1).data.numpy()
        preds = probility[:, 0] <= probility[:, 1]
        preds = preds.astype(int)[ind_sort]
        mols = np.asarray([i[-1].decode() for i in f['epoch_0000']['test']['mol']])[ind_sort]
        headers = np.asarray(['mol', 'predict_class', 'probility'])
        rank = np.asarray([i for i in range(1, len(mols)+1)])
        result = pd.DataFrame(np.vstack((rank, mols, preds, out)))
        result.index=['rank', 'pdb_name', 'predict_class', 'probility']
        result.to_csv(os.path.join(self.output_dir, 'deeprank.out'), encoding='utf-8', index=False)

#test the class
'''
pdb_source = '/trinity/login/xxu/data/test_CNN/1E6E_4'
outdir = '/trinity/login/xxu/data/test_CNN/1E6E_4/out'
CNN = CNN_score(pdb_source=pdb_source, chain1='A', chain2='B', output_dir=outdir)
#CNN.write_pssm()
#CNN.generate_pssm()
CNN.create_database()
CNN.test_CNN()
print(CNN.analysis_result())
'''
