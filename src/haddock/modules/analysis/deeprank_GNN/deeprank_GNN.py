
import os
import glob

from deeprank_gnn.GraphGenMP import GraphHDF5
from deeprank_gnn.NeuralNet import NeuralNet
from deeprank_gnn.ginet import GINet

import h5py
import torch.nn.functional as F

pretrained_model = '/trinity/login/xxu/software/haddock3/src/haddock/modules/analysis/deeprank_GNN/fold6_treg_yfnat_b128_e20_lr0.001_4.pt'

class GNN_score():
    def __init__(self, pdb_path, pssm_path, output_dir):
        self.pdb_path = pdb_path
        self.pssm_path = pssm_path
        self.output_dir = output_dir

    def generate_graph(self):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        GraphHDF5(pdb_path=self.pdb_path, pssm_path=self.pssm_path,
        graph_type='residue', outfile ='graph.hdf5', nproc=4)
    
    def predict(self):
        output_path = os.path.join(self.output_dir, 'prediction.hdf5')
        gnn = GINet
        database_test = glob.glob('*.hdf5')
        model = NeuralNet(database_test, gnn, pretrained_model = pretrained_model)
        model.test(threshold=None,hdf5=output_path)
        f = h5py.File(output_path)
        raw_outputs = f['epoch_0000']['test']['outputs'][()]
        predict_output = [ 1 if i >= 0.3 else 0 for i in raw_outputs]
        return predict_output


'''
#test the GNN_score class
pdb_path = '/trinity/login/xxu/software/Deeprank-GNN/example/data/pdb/1ATN'
pssm_path = '/trinity/login/xxu/software/Deeprank-GNN/example/data/pssm/1ATN'
output_dir = '/trinity/login/xxu/software/Deeprank-GNN/example/data/output'
GNN_score = GNN_score(pdb_path, pssm_path, output_dir)
GNN_score.generate_graph()
predict_output = GNN_score.predict()
print(predict_output)
'''