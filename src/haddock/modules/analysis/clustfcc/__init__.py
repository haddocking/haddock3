"""HADDOCK3 FCC clustering module"""
import logging
import os
from pathlib import Path

from fcc.scripts import calc_fcc_matrix, cluster_fcc

from haddock import FCC_path
from haddock.gear.config_reader import read_config
from haddock.libs.libparallel import Scheduler
from haddock.libs.libsubprocess import Job
from haddock.modules import BaseHaddockModule
from haddock.libs.libontology import Format, ModuleIO, PDBFile


logger = logging.getLogger(__name__)

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseHaddockModule):

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = False
        super().__init__(order, path, initial_params, cns_script)

    @classmethod
    def confirm_installation(cls):
        dcfg = read_config(DEFAULT_CONFIG)
        exec_path = Path(FCC_path, dcfg['executable'])

        if not os.access(exec_path, mode=os.F_OK):
            raise Exception(f'Required {str(exec_path)} file does not exist.')

        if not os.access(exec_path, mode=os.X_OK):
            raise Exception(f'Required {str(exec_path)} file is not executable')

        return

    def run(self, **params):
        logger.info("Running [clustfcc] module")

        super().run(params)

        contact_executable = Path(FCC_path, self.params['executable'])

        # Get the models generated in previous step
        models_to_cluster = [p for p in self.previous_io.output if p.file_type == Format.PDB]

        first_model = models_to_cluster[0]
        topologies = first_model.topology

        # Calculate the contacts for each model
        logger.info('Calculating contacts')
        contact_jobs = []
        for model in models_to_cluster:
            pdb_f = Path(model.path, model.file_name)
            contact_f = Path(self.path, model.file_name.replace('.pdb', '.con'))
            job = Job(
                pdb_f,
                contact_f,
                contact_executable,
                self.params['contact_distance_cutoff'],
                )
            contact_jobs.append(job)

        contact_engine = Scheduler(contact_jobs, ncores=self.params['ncores'])
        contact_engine.run()

        contact_file_l = []
        not_found = []
        for job in contact_jobs:
            if not job.output.exists():
                # NOTE: If there is no output, most likely the models are not in contact
                #  there is no way of knowing how many models are not in contact, it can be
                #  only one, or could be all of them.
                not_found.append(job.input.name)
                logger.warning(f'Contact was not calculated for {job.input.name}')
            else:
                contact_file_l.append(str(job.output))

        if not_found:
            # No contacts were calculated, we cannot cluster
            self.finish_with_error("Several files were not generated:"
                                   f" {not_found}")

        logger.info('Calculating the FCC matrix')
        parsed_contacts = calc_fcc_matrix.parse_contact_file(contact_file_l, False)

        # Imporant: matrix is a generator object, be careful with it
        matrix = calc_fcc_matrix.calculate_pairwise_matrix(parsed_contacts, False)

        # write the matrix to a file, so we can read it afterwards and don't
        #  need to reinvent the wheel handling this
        fcc_matrix_f = Path(self.path, 'fcc.matrix')
        with open(fcc_matrix_f, 'w') as fh:
            for data in list(matrix):
                data_str = f"{data[0]} {data[1]} {data[2]:.2f} {data[3]:.3f}"
                data_str += os.linesep
                fh.write(data_str)
        fh.close()

        # Cluster
        logger.info('Clustering...')
        pool = cluster_fcc.read_matrix(
            fcc_matrix_f,
            self.params['fraction_cutoff'],
            self.params['strictness'],
            )

        _, clusters = cluster_fcc.cluster_elements(
            pool,
            threshold=self.params['threshold'],
            )

        # Prepare output and read the elements
        cluster_dic = {}
        if clusters:
            # use fcc's output
            cluster_out = Path(self.path, 'cluster.out')
            with open(cluster_out, 'w') as fh:
                cluster_fcc.output_clusters(fh, clusters)
            fh.close()

            # Extract the cluster elements
            # Q: Can we do this without having to re-open the file?
            with open(cluster_out, 'r') as fh:
                for line in fh.readlines():
                    data = line.split()
                    cluster_id = int(data[1])
                    cluster_dic[cluster_id] = []
                    cluster_elements = list(map(int, data[4:]))
                    for element in cluster_elements:
                        pdb = models_to_cluster[element - 1]
                        cluster_dic[cluster_id].append(pdb)
        else:
            logger.warning('No clusters were found')

        # Save module information
        io = ModuleIO()
        io.add(models_to_cluster)
        io.add(cluster_dic, "o")
        io.save(self.path)
