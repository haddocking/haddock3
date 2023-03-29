import os
import h5py
import pandas as pd

from deeprank_gnn.GraphGenMP import GraphHDF5
from deeprank_gnn.NeuralNet import NeuralNet
from deeprank_gnn.ginet import GINet


class GNN_score:
    def __init__(self, pdb_path, output_dir):
        self.pdb_path = pdb_path
        self.output_dir = output_dir
        graph_path = self.generate_graph()
        output_path = self.predict(graph_path)
        self.analysis(output_path)

    def generate_graph(self):
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        graph_path = os.path.join(self.output_dir, "graph.hdf5")
        GraphHDF5(
            pdb_path=self.pdb_path,
            pssm_path=None,
            graph_type="residue",
            outfile=graph_path,
            nproc=10,
        )
        return graph_path

    def predict(self, graph_path):
        # graph_path = self.generate_graph()
        output_path = os.path.join(self.output_dir, "prediction.hdf5")
        gnn = GINet
        node_feature = ["type", "polarity", "bsa", "charge"]
        edge_feature = ["distance"]
        target = "fnat"
        pretrained_model = "/trinity/login/xxu/software/haddock3/src/haddock/modules/analysis/deeprank_GNN/treg_yfnat_b128_e10_lr0.001_8.pth.tar"
        #need to change the path
        model = NeuralNet(
            graph_path,
            gnn,
            edge_feature=edge_feature,
            node_feature=node_feature,
            target=target,
            pretrained_model=pretrained_model,
        )
        model.test(hdf5=output_path)
        return output_path

    def analysis(self, output_path):
        f = h5py.File(output_path, "r")
        csv_path = os.path.join(self.output_dir, "prediction.csv")
        mol_name = f["epoch_0000"]["test"]["mol"][()].astype(str)
        raw_outputs = f["epoch_0000"]["test"]["raw_outputs"][()]
        df = pd.DataFrame({"pdb": mol_name, "fnat": raw_outputs})
        sorted_df = df.sort_values(by="fnat", ascending=False)
        sorted_df.to_csv(csv_path, index=False)
