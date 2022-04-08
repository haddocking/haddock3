"""CapriJob file."""
from haddock import log


class CapriJob:
    """A Job dedicated to the caprieval module."""

    def __init__(
            self,
            output,
            params,
            capri_obj):

        core = capri_obj.core
        log.info(f"core {core}, initialising CAPRI...")
        log.info(f"core {core}, # of models : {len(capri_obj.model_list)}")
        self.output = output
        self.params = params
        self.capri_obj = capri_obj
        log.info(f"core {core}, CAPRI initialised")

    def run(self):
        """Run this CapriJob."""
        # adding models
        for struct in self.capri_obj.model_list:
            _ = self.capri_obj.add_chain_from_segid(struct.rel_path)
        
        # running calculations
        if self.params["fnat"]:
            log.debug(f"core {self.capri_obj.core}, calculating FNAT")
            fnat_cutoff = self.params["fnat_cutoff"]
            log.debug(f" cutoff: {fnat_cutoff}A")
            self.capri_obj.fnat(cutoff=fnat_cutoff)

        if self.params["irmsd"]:
            log.debug(f"core {self.capri_obj.core}, calculating I-RMSD")
            irmsd_cutoff = self.params["irmsd_cutoff"]
            log.debug(f" cutoff: {irmsd_cutoff}A")
            self.capri_obj.irmsd(cutoff=irmsd_cutoff)

        if self.params["lrmsd"]:
            log.debug(f"core {self.capri_obj.core}, calculating L-RMSD")
            self.capri_obj.lrmsd()

        if self.params["ilrmsd"]:
            log.debug(f"core {self.capri_obj.core}, calculating I-L-RMSD")
            ilrmsd_cutoff = self.params["irmsd_cutoff"]
            log.debug(f" cutoff: {ilrmsd_cutoff}A")

            self.capri_obj.ilrmsd(
                cutoff=ilrmsd_cutoff,
                )

        if self.params["dockq"]:
            log.debug(f"core {self.capri_obj.core}, calculating DockQ metric")
            self.capri_obj.dockq()
        
        self.capri_obj.output(
            self.params["clt_threshold"],
            sortby_key=self.params["sortby"],
            sort_ascending=self.params["sort_ascending"],
            )
        return
