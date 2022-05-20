"""Create Markdown pages for module's default parameters."""
import os
from collections.abc import Mapping
from pathlib import Path

from haddock import haddock3_repository_path
from haddock.modules import modules_category
from haddock.modules.util import read_all_default_configs_yaml


EXPLEVELS = (
    'easy',
    'expert',
    'guru',
    'hidden',
    )


# prepare YAML markdown files
def prepare_yaml_markdown():
    """
    Prepare Markdown files to document YAML default files in modules.

    These files are written to the docs/ folder but not stagged to github.
    They are used only by Sphinx to generate the HTML documentation pages.
    """
    def do_text(name, param):
        """Create text from parameter dictionary."""
        text = [
            f'{name}',
            f'{"`" * len(name)}',
            '',
            f'| *default*: {param["default"]!r}',
            f'| *type*: {param["type"]}',
            f'| *title*: {param["title"]}',
            f'| *short*: {param["short"]}',
            f'| *long*: {param["long"]}',
            f'| *group*: {param["group"]}',
            f'| *explevel*: {param["explevel"]}',
            '',
            ]

        return os.linesep.join(text)

    def loop_params(config):
        """Treat parameters for module."""
        for name, data in config.items():
            if isinstance(data, Mapping) and "default" not in data:
                loop_params(data)
                continue

            elif isinstance(data, Mapping):
                explevel = data["explevel"]

                text = do_text(name, data)

                if explevel == 'easy':
                    easy.append(text)
                elif explevel == 'expert':
                    expert.append(text)
                elif explevel == 'guru':
                    guru.append(text)
                elif explevel == 'hidden':
                    continue
                else:
                    emsg = f'explevel {explevel!r} is not expected'
                    raise AssertionError(emsg)

        easy.append('')
        expert.append('')
        guru.append('')

    configs = read_all_default_configs_yaml()

    for module, params in configs.items():
        easy = ['Easy', '----', '']
        expert = ["Expert", '------', '']
        guru = ['Guru', '----', '']

        loop_params(params)

        category = modules_category[module]
        params_folder = Path(
            haddock3_repository_path,
            'docs',
            'modules',
            category,
            'params',
            )
        params_folder.mkdir(exist_ok=True)

        doc = []
        for list_ in (easy, expert, guru):
            if len(list_) > 4:
                doc.extend(list_)

        with open(Path(params_folder, f'{module}.rst'), 'w') as fout:
            fout.write(os.linesep + os.linesep)
            fout.write(os.linesep.join(doc))


prepare_yaml_markdown()
