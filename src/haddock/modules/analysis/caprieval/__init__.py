"""Calculate CAPRI metrics."""
from pathlib import Path

from haddock.libs.libontology import Format, ModuleIO
from haddock.modules import BaseHaddockModule
from haddock.modules.analysis.caprieval.capri import CAPRI


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to calculate the CAPRI metrics."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, *ignore, init_params=DEFAULT_CONFIG,
                 **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm if contact executable is compiled."""
        return

    def _run(self):
        """Execute module."""
        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            _e = "This module cannot come after one that produced an iterable."
            self.finish_with_error(_e)

        models_to_calc = [
            p for p in self.previous_io.output if p.file_type == Format.PDB
            ]

        #  by default modes_to_calc should have been sorted by the module
        #  that produced it
        best_model = Path(models_to_calc[0].path, models_to_calc[0].file_name)

        if self.params["reference"]:
            reference = Path(self.params["reference"])
        else:
            self.log(
                "No reference was given. "
                "Using best ranking structure from previous step")
            reference = best_model

        try:
            capri = CAPRI(
                reference,
                models_to_calc,
                receptor_chain=self.params["receptor_chain"],
                ligand_chain=self.params["ligand_chain"],
                aln_method=self.params["alignment_method"],
                path=self.path
                )
        except Exception as e:
            self.finish_with_error(e)

        if self.params["fnat"]:
            self.log("Calculating FNAT")
            fnat_cutoff = self.params["fnat_cutoff"]
            self.log(f" cutoff: {fnat_cutoff}A")
            capri.fnat(cutoff=fnat_cutoff)

        if self.params["irmsd"]:
            self.log("Calculating I-RMSD")
            irmsd_cutoff = self.params["irmsd_cutoff"]
            self.log(f" cutoff: {irmsd_cutoff}A")
            capri.irmsd(cutoff=irmsd_cutoff)

        if self.params["lrmsd"]:
            self.log("Calculating L-RMSD")
            capri.lrmsd()

        if self.params["ilrmsd"]:
            self.log("Calculating I-L-RMSD")
            ilrmsd_cutoff = self.params["irmsd_cutoff"]
            self.log(f" cutoff: {ilrmsd_cutoff}A")

            capri.ilrmsd(
                cutoff=ilrmsd_cutoff,
                )

        output_fname = Path(self.path, "capri.tsv")
        self.log(f" Saving output to {output_fname.name}")
        capri.output(
            output_fname,
            sortby_key=self.params["sortby"],
            sort_ascending=self.params["sort_ascending"],
            rankby_key=self.params["rankby"],
            rank_ascending=self.params["sort_ascending"],
            )

        selected_models = models_to_calc
        io = ModuleIO()
        io.add(selected_models, "o")
        io.save(self.path)
