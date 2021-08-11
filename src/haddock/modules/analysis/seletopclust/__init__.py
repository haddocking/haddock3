"""HADDOCK3 module to select a top cluster/model"""
import logging
from pathlib import Path
from haddock.modules import BaseHaddockModule
from haddock.ontology import ModuleIO

logger = logging.getLogger(__name__)

class HaddockModule(BaseHaddockModule):

    def __init__(self, order, path, *ignore, **everything):
        recipe_path = Path(__file__).resolve().parent
        cns_script = ''
        defaults = recipe_path / "seletopclust.toml"
        super().__init__(order, path, cns_script, defaults)

    def run(self, **params):
        logger.info("Running [selecttopclust] module")

        # Pool of jobs to be executed by the CNS engine
        jobs = []

        # Get the models generated in previous step
        if not type(self.previous_io) == iter:
            # this module needs to come after one that produced an iterable
            pass

        # retrieve the clusters from a dictionary generated in the previous step
        #  the cluster_id just tells us how populated a given cluster is
        # if we discard this value we can have lists
        average_score_dic = {}
        cluster_dic = {}
        # Q: Why this [0] here?
        for cluster_id in self.previous_io.output[0]:
            cluster_id = int(cluster_id)
            # sort the models inside the cluster based on its score
            # TODO: refactor this, its ugly :p
            list_to_be_sorted = [(e, e.score) for e in self.previous_io.output[0][str(cluster_id)]]
            list_to_be_sorted.sort(key = lambda x: x[1])
            structure_list = [e[0] for e in list_to_be_sorted]
            cluster_dic[cluster_id] = structure_list

            # get the average score of the cluster based on ALL the elements
            scores = [e[1] for e in list_to_be_sorted]
            average_score = sum(scores) / float(len(scores))
            average_score_dic[cluster_id] = average_score
            
        # sort the clusters based on their average
        sorted_dic = dict(sorted(average_score_dic.items(), key=lambda item: item[1]))

        # =====
        # TODO: implement behaviour if top_cluster is a list
        #  this way we could for example only run flexible refinement for the
        #  top4 clusters
        # =====

        # which cluster should we retrieve?
        # top_cluster = 1 == the best one, should be index 0
        target_cluster_id = list(sorted_dic.keys())[params['top_cluster'] - 1]

        # how many models should we output?
        if params['top_models'] == 'all':
            models_to_output = cluster_dic[target_cluster_id]
        else:
            models_to_output = cluster_dic[target_cluster_id][:params['top_models']]

        # Save module information
        io = ModuleIO()
        io.add(models_to_output, "o")
        io.save(self.path)