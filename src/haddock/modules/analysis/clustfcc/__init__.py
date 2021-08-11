"""HADDOCK3 FCC clustering module"""
import logging
from os import linesep
from pathlib import Path
import multiprocessing
from haddock.modules import BaseHaddockModule
from haddock.engine import Job, Engine
from haddock.cns.util import generate_default_header, load_ambig
from haddock.cns.util import load_workflow_params, prepare_multiple_input
from haddock.ontology import Format, ModuleIO, PDBFile
from fcc.scripts import calc_fcc_matrix, make_contacts, cluster_fcc

logger = logging.getLogger(__name__)

# def calc_contact(executable, pdb_f, cutoff=5.0):
#     # fcc-wrapper
#     contact_f = Path(pdb_f.replace('.pdb', '.contacts'))
#     make_contacts._calculate_contacts(executable, pdb_f, str(cutoff))
#     if contact_f.exists():
#     # if os.path.isfile(contact_f):
#         return contact_f
#     else:
#         return ''


# def cluster_fcc(structures, contact_executable, cutoff=0.75, np=1):
    # """Use FCC to cluster structures."""
    # logger.info()

    #     # Calculate contacts
    #     ga_log.info('FCC - Calculating contacts')
    # pool = multiprocessing.Pool(np)    #     input_structure_l = self.structure_list[:top]
    #     for pdb in input_structure_l:
    #         pool.apply_async(self.calc_contact, args=(self.contact_executable,
    #                                                   pdb))

    #     pool.close()
    #     pool.join()

    #     contact_file_l = []
    #     for pdb in input_structure_l:
    #         contact_f = pdb.replace('.pdb', '.contacts')
    #         if os.path.isfile(contact_f):
    #             contact_file_l.append(contact_f)

    #     if not contact_file_l:
    #         ga_log.warning('No contacts were calculated')

    #     # Calculate matrix
    #     ga_log.info('FCC - Calculating matrix')
    #     parsed_contacts = calc_fcc_matrix.parse_contact_file(contact_file_l,
    #                                                          False)

    #     # matrix is a generator object, be careful with it
    #     matrix = calc_fcc_matrix.calculate_pairwise_matrix(parsed_contacts,
    #                                                        False)

    #     # write it to a file, so we can read it afterwards and don't
    #     #  need to reinvent the wheel
    #     fcc_matrix_f = f'{self.analysis_path}/fcc.matrix'
    #     with open(fcc_matrix_f, 'w') as fh:
    #         for data in list(matrix):
    #             data_str = f"{data[0]} {data[1]} {data[2]:.2f} {data[3]:.3f}"
    #             data_str += os.linesep
    #             fh.write(data_str)
    #     fh.close()

    #     # cluster
    #     ga_log.info('FCC - Clustering')
    #     pool = cluster_fcc.read_matrix(fcc_matrix_f, cutoff, strictness=0.75)
    #     _, clusters = cluster_fcc.cluster_elements(pool, 4)

    #     if clusters:
    #         ga_log.info(f'FCC - {len(clusters)} clusters identified')
    #         # use fcc's output
    #         cluster_out = f'{self.analysis_path}/cluster.out'
    #         with open(cluster_out, 'w') as fh:
    #             cluster_fcc.output_clusters(fh, clusters)
    #         fh.close()

    #         # read it again!
    #         with open(cluster_out, 'r') as fh:
    #             for line in fh.readlines():
    #                 data = line.split()
    #                 cluster_id = int(data[1])
    #                 self.cluster_dic[cluster_id] = []
    #                 cluster_elements = list(map(int, data[4:]))
    #                 for element in cluster_elements:
    #                     structure_name = self.structure_list[element - 1]
    #                     element_name = pathlib.Path(structure_name).stem
    #                     self.cluster_dic[cluster_id].append(element_name)
    #     else:
    #         ga_log.info('FCC - No clusters identified')

    #     return self.cluster_dic

    # Calculate contacts

    # Calculate Matrix

    # Cluster
    # pass


class HaddockModule(BaseHaddockModule):

    def __init__(self, order, path, *ignore, **everything):
        recipe_path = Path(__file__).resolve().parent
        cns_script = False
        defaults = recipe_path / "clustfcc.toml"
        super().__init__(order, path, cns_script, defaults)

    def run(self, **params):
        logger.info("Running [clustfcc] module")

        # Q: Here we need to find haddock3/src/fcc, is there a better way of doign it?
        contact_executable = Path(self.defaults_path.parent.parent.parent.parent.parent.parent,
                                  self.defaults['params']['executable'])

        # # Pool of jobs to be executed by the CNS engine
        # jobs = []

        # Get the models generated in previous step
        models_to_cluster = [p for p in self.previous_io.output if p.file_type == Format.PDB]

        first_model = models_to_cluster[0]
        topologies = first_model.topology

        contact_jobs = []
        contact_cutoff = params['cutoff']
        for model in models_to_cluster:
            pdb_f = Path(model.path, model.file_name)
            contact_f = Path(self.path, model.file_name.replace('.pdb', '.con'))
            job = Job(pdb_f, contact_f, contact_executable, contact_cutoff)
            contact_jobs.append(job)

        contact_engine = Engine(contact_jobs)
        contact_engine.run()
        print('idk')



            #     input_structure_l = self.structure_list[:top]
    #     for pdb in input_structure_l:
    #         pool.apply_async(self.calc_contact, args=(self.contact_executable,
    #                                                   pdb))

    #     pool.close()
    #     pool.join()

    #     contact_file_l = []
    #     for pdb in input_structure_l:
    #         contact_f = pdb.replace('.pdb', '.contacts')
    #         if os.path.isfile(contact_f):
    #             contact_file_l.append(contact_f)

    #     if not contact_file_l:
    #         ga_log.warning('No contacts were calculated')


        cluster_dic = cluster_fcc(structures=models_to_cluster,
                                  contact_executable=contact_executable,
                                  cutoff=params['cutoff'])

        # structure_list = []
        # for idx in range(params['sampling']):
        #     inp_file = generate_docking(
        #         idx,
        #         models_to_dock,
        #         self.path,
        #         self.recipe_str,
        #         self.defaults,
        #         ambig=params.get('ambig', None),
        #         )

        #     out_file = self.path / f"rigidbody_{idx}.out"
        #     structure_file = self.path / f"rigidbody_{idx}.pdb"
        #     structure_list.append(structure_file)

        #     job = CNSJob(inp_file, out_file, cns_folder=self.cns_folder_path)

        #     jobs.append(job)

        # # Run CNS engine
        # logger.info(f"Running CNS engine with {len(jobs)} jobs")
        # engine = CNSEngine(jobs)
        # engine.run()
        # logger.info("CNS engine has finished")

        # Check for generated output, fail it not all expected files are found
        expected = []
        not_found = []
        for model in structure_list:
            if not model.exists():
                not_found.append(model.name)
            pdb = PDBFile(model, path=self.path)
            pdb.topology = topologies
            expected.append(pdb)
        if not_found:
            self.finish_with_error("Several files were not generated:"
                                   f" {not_found}")

        # Save module information
        io = ModuleIO()
        io.add(structure_list)
        io.add(expected, "o")
        io.save(self.path)