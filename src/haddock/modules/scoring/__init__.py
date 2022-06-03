"""HADDOCK3 modules to score modules."""
from os import linesep

from haddock.modules.base_cns_module import BaseCNSModule


class ScoringModule(BaseCNSModule):
    """Parent class for Scoring modules."""

    def output(self, output_fname, sep="\t"):
        """Save the output in comprehensive tables."""
        # prepares header
        header = sep.join(
            ("structure", "original_name", "md5", "score")
            ) + linesep

        # prepares a text generator
        text_generator = (
            f"{pdb.file_name}\t{pdb.ori_name}\t{pdb.md5}\t{pdb.score}"
            for pdb in self.output_models
            )

        # writes to disk only once
        with open(output_fname, "w") as fh:
            fh.write(header + linesep.join(text_generator))

        return
