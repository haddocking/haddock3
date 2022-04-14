"""Traces the input PDBs accross the modules."""
import os


class Traceback:
    """Class that represents the PDB traceback."""

    def __init__(self, io_data):
        self.model_data, self.module_list = self.read_io_data(io_data)

    @staticmethod
    def read_io_data(io_data):
        """Read the IO data from modules and organize it per-model."""
        model_dict = {}
        module_list = []
        for module_name in io_data:
            if '0_topoaa' in module_name:
                # ignore topoaa
                continue
            module_list.append(module_name)
            model_container = io_data[module_name].retrieve_models(
                individualize=True)
            for model in model_container:
                if model.ori_name not in model_dict:
                    model_dict[model.ori_name] = {}
                model_dict[model.ori_name][module_name] = model.score

        return model_dict, module_list

    def dump(self, output_f=None):
        """Dump the traceback info to a file."""
        header = 'model_name\t' + "\t".join(self.module_list)
        with open(output_f, 'w') as fh:
            fh.write(header + os.linesep)
            for model in self.model_data:
                row = ""
                row += f"{model}\t"
                for module in self.module_list:
                    try:
                        value = self.model_data[model][module]
                        value_str = f"{value:.2f}"
                    except KeyError:
                        value_str = '-'
                    row += f"{value_str}\t"
                fh.write(row + os.linesep)
