"""HADDOCK3 module to select a top cluster/model"""
import logging
from pathlib import Path
from haddock.modules import BaseHaddockModule
from haddock.libs.libontology import ModuleIO

logger = logging.getLogger(__name__)

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseHaddockModule):

    def __init__(
            self,
            order,
            path,
            *ignore,
            init_params=DEFAULT_CONFIG,
            **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        return

    def run(self, **params):
        logger.info("Running [seletopclusts] module")

        super().run(params)

        # Get the models generated in previous step
        if not type(self.previous_io) == iter:
            # this module needs to come after one that produced an iterable
            pass

        # retrieve the clusters from a dictionary generated in the previous
        # step the cluster_id just tells us how populated a given cluster is
        # if we discard this value we can have lists
        average_dic = {}
        cluster_dic = {}
        # Q: Why this [0] here?
        for cluster_id in self.previous_io.output[0]:
            cluster_id = int(cluster_id)
            # sort the models inside the cluster based on its score
            # TODO: refactor this, its ugly :p
            list_to_be_sorted = [(e, e.score) for e in self.previous_io.output[0][str(cluster_id)]]
            list_to_be_sorted.sort(key=lambda x: x[1])
            structure_list = [e[0] for e in list_to_be_sorted]
            cluster_dic[cluster_id] = structure_list

            # get the average score of the cluster based on ALL the elements
            scores = [e[1] for e in list_to_be_sorted]
            average_score = sum(scores) / float(len(scores))
            average_dic[cluster_id] = average_score

        # sort the clusters based on their average
        sorted_dic = sorted(average_dic.items(), key=lambda item: item[1])
        sorted_dic = dict(sorted_dic)

        # how many models should we output?
        models = []
        for select_id in self.params['top_cluster']:
            # which cluster should we retrieve?
            # top_cluster = 1 == the best one, should be index 0
            try:
                target_id = list(sorted_dic.keys())[select_id - 1]
            except IndexError:
                logger.warning(f'Cluster ranking #{select_id} not found,'
                               ' skipping selection')
                continue

            if self.params['top_models'] == 'all':
                for pdb in cluster_dic[target_id]:
                    models.append(pdb)
            else:
                for pdb in cluster_dic[target_id][:self.params['top_models']]:
                    models.append(pdb)

        # Save module information
        io = ModuleIO()
        io.add(models, "o")
        io.save(self.path)
