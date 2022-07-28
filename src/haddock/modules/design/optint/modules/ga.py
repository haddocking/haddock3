"""Genetic Algorithm module."""
import multiprocessing
import os
import random
import sys
import tempfile
from itertools import repeat
from operator import attrgetter
from pathlib import Path

import numpy as np
from deap import base, creator, tools
from pdbtools.pdb_tidy import tidy_pdbfile

from haddock import log
from haddock.modules.design.optint.modules.fitness import calc_haddockscore
from haddock.modules.design.optint.modules.pdb import (
    get_interface_dict,
    pdb_to_dict,
    )
from haddock.modules.design.optint.modules.residues import (
    CODE_TO_AA,
    THREE_TO_ONE,
    )


try:
    from collections.abc import Sequence
except ImportError:
    from collections import Sequence

creator.create("FitnessSingle", base.Fitness, weights=(-1.0,))
creator.create("Individual", list, fitness=creator.FitnessSingle)


ATOMS_TO_BE_MUTATED = ["C", "N", "CA", "O", "CB"]


class GA:
    """Genetic Algorithm Class."""

    def __init__(self, identifier, input_pdb, params):
        self.pdb = input_pdb
        self.id = identifier
        self.target_chain = params["target_chain"]
        self.random_seed = params["random_seed"]
        self.max_ngen = params["number_of_generations"]
        self.popsize = params["population_size"]
        self.cxpb = params["cxpb"]
        self.mutpb = params["mutpb"]
        self.eta = params["eta"]
        self.indpb = params["indpb"]
        self.number_of_models = params["output_models"]
        self.region = params["region"]
        self.initial_mutation = params["initial_mutation"]
        self.toolbox = None
        self.generation_dic = {}
        self.region_dict = self.get_region_to_be_optimized()
        if params["ncores"] > self.popsize:
            self.nproc = self.popsize
        else:
            self.nproc = params["ncores"]

    def setup(self):
        """Run the toolbox setup."""
        toolbox = base.Toolbox()
        toolbox.register(
            "attr",
            self.generate_individual,
            region_dict=self.region_dict,
            target_chain=self.target_chain,
            initial_mutation=self.initial_mutation,
            )
        toolbox.register(
            "individual", tools.initIterate, creator.Individual, toolbox.attr
            )
        toolbox.register("population", tools.initRepeat,
                         list, toolbox.individual)

        toolbox.register("mate", tools.cxTwoPoint)

        toolbox.register(
            "mutate_aa",
            self.custom_mutUniformInt,
            low=1,
            up=20,
            indpb=self.indpb,
            )

        toolbox.register("select", self.custom_selTournament, tournsize=3)

        toolbox.register("evaluate", self.fitness_function, self.pdb)

        self.toolbox = toolbox

    def run(self):
        """Run the genetic algorithm."""
        self.setup()
        pool = multiprocessing.Pool(processes=self.nproc)
        self.toolbox.register("map", pool.map)

        log.info("Genetic Algorithm parameters: ")
        log.info(f"  + population_size: {self.popsize}")
        log.info(f"  + max_number_of_generations: {self.max_ngen}")
        log.info(f"  + crossover_probability: {self.cxpb}")
        log.info(f"  + mutation_probability: {self.mutpb}")
        log.info(f"  + eta (crowding degree of mut): {self.eta}")
        log.info(f"  + indpb (independent prob): {self.indpb}")

        random.seed(self.random_seed)
        run = True
        log.info(f"Starting the simulation with {self.nproc} processors")
        log.info(f"  random_seed = {self.random_seed}")
        log.info("##########################################")
        log.info("Gen ... mean +-  sd   (min, max) unique_ind")
        pop = self.toolbox.population(n=self.popsize)
        ngen = 1
        while run:
            self.generation_dic[ngen] = {}

            offspring = self.toolbox.select(pop, len(pop))

            # Clone the population created
            offspring = list(map(self.toolbox.clone, offspring))

            # Add generation/individual information
            if any([type(e) == tuple for e in offspring]):
                offspring = [(ngen, i, e[2], self.id, self.target_chain)
                             for i, e in enumerate(offspring)]
            else:
                offspring = [(ngen, i, e, self.id, self.target_chain)
                             for i, e in enumerate(offspring)]

            # Apply crossover on the offspring
            for element_1, element_2 in zip(offspring[::2], offspring[1::2]):
                child1 = element_1[2]
                child2 = element_2[2]
                if random.random() < self.cxpb:  # nosec
                    self.toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            # Apply mutation on the offspring
            for element in offspring:
                # mutant = element[2]
                if random.random() < self.mutpb:  # nosec
                    self.toolbox.mutate_aa(element)
                    del element[2].fitness.values

            # evaluate only the mutated and crossed
            invalid_ind = [ind for ind in offspring if not ind[2].fitness.valid]
            try:
                fitnesses = self.toolbox.map(self.toolbox.evaluate, invalid_ind)
            except Exception as e:
                log.warning(e)
                log.error("Fitness function could not be executed")
                sys.exit()

            for ind, fit in zip(invalid_ind, fitnesses):
                ind[2].fitness.values = fit  # , fit

            # replace the old population by the offspring
            pop[:] = offspring

            fitness_values = []
            for idx, ind in enumerate(pop):
                fitness_v = ind[2].fitness.values
                self.generation_dic[ngen][idx] = {
                    "individual": ind,
                    "fitness": fitness_v,
                    }

                fitness_values.append(ind[2].fitness.values[0])

            fitness_summary = self.summary(fitness_values)

            ngen_str = str(ngen).rjust(3, "0")
            log.info(
                f"Gen {ngen_str} {fitness_summary['mean']:.2f} "
                f"+- {fitness_summary['std']:.2f} "
                f"({fitness_summary['min']:.2f}, "
                f"{fitness_summary['max']:.2f}) "
                f"eval={len(invalid_ind)}"
                )

            if ngen == self.max_ngen:
                log.info(
                    f"Simulation reached maximum number of "
                    f"generations, stopping at {ngen}."
                    )
                run = False

            ngen += 1

        pool.close()
        pool.join()

        # ga complete, get the best models
        best_pdb_l = self.retrieve_models()

        return best_pdb_l

    def retrieve_models(self):
        """Get the best PDB models."""
        results = []

        for f in Path(".").glob("*.score"):
            with open(f) as fh:
                data = fh.readlines()[0].split()
                pdb = data[0]
                psf = pdb.replace(".pdb", ".psf")
                score = float(data[3])
                results.append((pdb, psf, score))

        results.sort(key=lambda x: x[2])

        return results[: self.number_of_models]

    def get_region_to_be_optimized(self):
        """Get the region to be optimized."""
        if self.region:
            pdb_dict = pdb_to_dict(str(self.pdb.rel_path))
            for resnum in pdb_dict[self.target_chain]:
                if resnum not in self.region:
                    pdb_dict[self.target_chain][resnum] = None
        else:
            log.info("Getting interfacial residues...")
            pdb_dict = get_interface_dict(str(self.pdb.rel_path))

        return pdb_dict

    @staticmethod
    def fitness_function(pdb_f, element):
        """Fitness function."""
        ngen, ident, individual, identifier, target_chain = element
        cwd = Path.cwd()

        # 1 - Use the chromossome to generate the mutated PDB
        mut_pdb_l = []
        pos = -1
        current_resnum = []
        with open(pdb_f.rel_path, "r") as fh:
            for line in fh.readlines():
                if line.startswith("ATOM"):
                    chain = line[21]
                    resnum = int(line[22:26])
                    atom_name = line[12:16].strip()
                    if target_chain == chain:
                        if resnum not in current_resnum:
                            current_resnum.append(resnum)
                        pos = len(current_resnum) - 1
                        target_pos = individual[pos]
                        if target_pos is not None:
                            mut_resname = CODE_TO_AA[target_pos]
                            if atom_name in ATOMS_TO_BE_MUTATED:
                                # mutate
                                mut_line = line[:17] + \
                                    mut_resname + line[20:]
                                mut_pdb_l.append(mut_line)
                        else:
                            # ignore this atom since its a side chain
                            mut_pdb_l.append(line)
                    if target_chain != chain:
                        mut_pdb_l.append(line)

        input_mutated_pdb = tempfile.NamedTemporaryFile(
            delete=False, mode="w", dir=str(cwd), suffix=".pdb"
            )
        input_mutated_pdb.write("".join(mut_pdb_l))
        input_mutated_pdb.close()
        mutated_pdb_fname = Path(input_mutated_pdb.name)

        # 2 - Tidy up the PDB
        tidy_pdb_l = []
        with open(mutated_pdb_fname, "r") as inp_fh:
            for line in tidy_pdbfile(inp_fh):
                tidy_pdb_l.append(line)

        tidy_mutated_pdb = tempfile.NamedTemporaryFile(
            delete=False, mode="w", dir=str(cwd), suffix=".pdb"
            )

        tidy_mutated_pdb.write("".join(tidy_pdb_l))
        tidy_mutated_pdb.close()
        tidy_model_path = Path(tidy_mutated_pdb.name)

        # 3 - Calculate Fitness
        fitness, output_pdb = calc_haddockscore(tidy_model_path)

        # 4 - Dump information for later analysis
        ngen_str = str(ngen).zfill(3)
        ident_str = str(ident).zfill(3)

        individual_pdb_f = Path(
            f"optint_{identifier}_{ngen_str}_{ident_str}.pdb"
            )
        individual_psf_f = Path(
            f"optint_{identifier}_{ngen_str}_{ident_str}.psf"
            )

        output_psf = output_pdb.with_suffix('.psf')

        output_pdb = output_pdb.rename(individual_pdb_f)
        output_psf = output_psf.rename(individual_psf_f)

        sequence = "".join(
            [THREE_TO_ONE[CODE_TO_AA[j]] for j in individual if j is not None]
            )

        score_f = output_pdb.with_suffix(".score")
        with open(score_f, "w") as f:
            f.write(
                f"{output_pdb}\t{ngen}\t{identifier}\t{ident}\t"
                f"{fitness}\t{sequence}" + os.linesep
                )

        # 5 - Clean
        mutated_pdb_fname.unlink()
        tidy_model_path.unlink()

        # 6 - Return
        return [fitness]

    @staticmethod
    def generate_individual(region_dict, target_chain, initial_mutation):
        """Generate the individual."""
        individual = []
        for resnum in region_dict[target_chain]:
            value = region_dict[target_chain][resnum]
            if value is not None:
                if initial_mutation == "random":
                    target_mutation = random.randint(1, 20)
                elif initial_mutation == "alanine":
                    target_mutation = 1  # ALA
                elif initial_mutation == "keep":
                    target_mutation = value  # do nothing
                individual.append(target_mutation)
            else:
                individual.append(None)
        return individual

    @staticmethod
    def summary(value_list):
        """Calculate the summary statistics of a value list."""
        mean = np.mean(value_list)
        std = np.std(value_list)
        max_v = max(value_list)
        min_v = min(value_list)
        return {"mean": mean, "std": std, "max": max_v, "min": min_v}

    @staticmethod
    def custom_mutUniformInt(individual, low, up, indpb):
        """Generate mutated individuals."""
        # tools.mutUniformInt
        size = len(individual[2])
        if not isinstance(low, Sequence):
            low = repeat(low, size)
        elif len(low) < size:
            raise IndexError(
                "low must be at least the size of individual: %d < %d"
                % (len(low), size)
                )
        if not isinstance(up, Sequence):
            up = repeat(up, size)
        elif len(up) < size:
            raise IndexError(
                "up must be at least the size of individual: %d < %d" % (
                    len(up), size)
                )
        for i, xl, xu in zip(range(size), low, up):
            if individual[2][i] is not None:
                if random.random() < indpb:
                    individual[2][i] = random.randint(xl, xu)
        return (individual,)

    @staticmethod
    def custom_selTournament(individuals, k, tournsize, fit_attr="fitness"):
        """Select the best individuals."""
        # tools.selTournament
        chosen = []
        for _ in range(k):
            # Randomly select 3 individuals
            aspirants = tools.selRandom(individuals, tournsize)
            # Find out who is the best
            try:
                best_individual = max(aspirants, key=attrgetter(fit_attr))
            except AttributeError:
                # This individual has the generation and identification
                #  information
                ind_dict = dict((j, e[2].fitness)
                                for j, e in enumerate(aspirants))
                sorted_ind_tuple = sorted(
                    ind_dict.items(), key=lambda item: item[1])
                best_individual_key = sorted_ind_tuple[0][0]
                best_individual = aspirants[best_individual_key]

            chosen.append(best_individual)

        return chosen
