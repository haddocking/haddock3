"""
Generate the docs/src/ API-reference pages for HADDOCK3.

Run automatically by docs/conf.py (via Sphinx's `builder-inited` hook)
right after `sphinx-apidoc` has dropped its flat, auto-generated `.rst`
stubs into docs/src/. This script does two things with them:

1. Renders each module's `defaults.yaml` parameters into a `.rst` page
   (docs/src/modules/<category>/params/<module>.rst), and appends an
   `.. include::` for it to the module's own stub, so the parameters show
   up on that module's documentation page.
2. Moves and retitles the apidoc stubs into the docs/src/modules/,
   docs/src/reference/, and docs/src/clients/ trees, using the friendly
   titles from titles.yaml (next to this script) instead of the raw
   dotted module names autodoc would otherwise use as page titles.

Everything under docs/src/ is generated output: it is gitignored and
rebuilt from scratch on every docs build, never committed to the repo.
"""

import os
import shutil
from collections.abc import Mapping
from pathlib import Path

from haddock.core.typing import ParamMap
from haddock.libs.libio import read_from_yaml

HADDOCK_REL_SRC_PATH = Path("src", "haddock")
DOCS_SRC = Path("docs", "src")

# Human-readable page titles live in titles.yaml, not here -- see that file
# for how to add a title for a new module/category/reference/CLI entry.
_TITLES = read_from_yaml(Path(__file__).parent / "titles.yaml")
MODULE_TITLE_DICT = _TITLES["modules"]
CATEGORY_TITLE_DICT = _TITLES["categories"]
REFERENCE_TITLE_DICT = _TITLES["reference"]
CLI_TITLE_DICT = _TITLES["clis"]


class HeadingController:
    """
    Control headings.

    reStructured text headings are defined by punctuation characters.

    In HADDOCK3 docs we use the order: '=', '-', '`', '~', '*'.

    The first heading tags is taken by the main docs. Therefore,
    `HeadingController` manages only from the second ('-') onward.

    Read more at: https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html#headings
    """  # noqa: E501

    def __init__(self) -> None:
        self.title_headings = ["-", "`", "~", "*"]
        self._idx = 0

    @property
    def next(self) -> str:
        """Give the next heading char."""
        return self.title_headings[self._idx + 1]

    @property
    def current(self) -> str:
        """Give the current heading char."""
        return self.title_headings[self._idx]

    def reset(self) -> None:
        """Reset to the first heading."""
        self._idx = 0

    def increase(self) -> None:
        """Increase current heading."""
        self._idx += 1


HEADING = HeadingController()


def change_title(rst_file: Path, title: str) -> None:
    """
    Change the title of the rst file.

    Parameters
    ----------
    rst_file : Path
        Path to the rst file.
    title : str
        New title.
    """
    with open(rst_file, "r") as fin:
        lines = fin.readlines()
    with open(rst_file, "w") as fout:
        for ln, line in enumerate(lines):
            if ln == 0:
                line = title + os.linesep
            else:
                line = line
            fout.write(line)


def process_category_file(category: str) -> None:
    """
    Process the category file.

    Parameters
    ----------
    category : str
        Category name.
    """
    category_rst = Path(DOCS_SRC, f"haddock.modules.{category}.rst")

    target_path = Path(DOCS_SRC, "modules", category)
    # make category folder if it does not exist
    target_path.mkdir(exist_ok=True, parents=True)
    target_rst = Path(target_path, "index.rst")
    shutil.move(category_rst, target_rst)
    # change title
    if category in CATEGORY_TITLE_DICT:
        title = CATEGORY_TITLE_DICT[category]
        change_title(target_rst, title)


def process_module_file(category: str, module_name: str) -> None:
    """
    Process the module file.

    Parameters
    ----------
    category : str
        Category name.
    module_name : str
        Module name.
    """
    module_rst = Path(DOCS_SRC, f"haddock.modules.{category}.{module_name}.rst")
    target_rst = Path(DOCS_SRC, "modules", category, module_rst.name)
    shutil.move(module_rst, target_rst)
    # does the submodule exist?
    submodule_gen = DOCS_SRC.glob(f"haddock.modules.{category}.{module_name}.*.rst")
    submodule_list = list(submodule_gen)
    if len(submodule_list) != 0:
        submodule_rst = submodule_list[0]
        submodule_target_rst = Path(DOCS_SRC, "modules", category, submodule_rst.name)
        shutil.move(submodule_rst, submodule_target_rst)

    with open(target_rst, "a") as fout:
        fout.write(
            f"{os.linesep}Default Parameters{os.linesep}"
            f"---------------{os.linesep}"
            f".. include:: params/{module_name}.rst" + os.linesep + os.linesep
        )
    # change title
    if module_name in MODULE_TITLE_DICT:
        title = MODULE_TITLE_DICT[module_name]
        change_title(target_rst, title)


# prepare YAML markdown files
def main() -> None:
    """
    Prepare restructured text files from YAML default configs in modules.

    These files are written to the 'docs/src/' folder but not stagged to
    github. Instead, they are used only by Sphinx to generate the HTML
    documentation pages.
    """
    # uses this pattern instead of importing:
    # from haddock.modules import modules_category
    # to avoid importing dependencies of the haddock modules packages
    pattern = Path("modules", "*", "*", "*.yaml")
    configs = HADDOCK_REL_SRC_PATH.glob(str(pattern))
    processed_categories = []
    # create RST pages for all modules' configuration files.
    for config in configs:
        if "_template" in str(config):
            continue

        module_name = config.parents[0].name
        category = config.parents[1].name
        params = read_from_yaml(config)

        # ignore empty modules - currently topocg for example
        if len(params) == 0:
            continue
        # if the category has not been processed yet, copy the category file
        if category not in processed_categories:
            process_category_file(category)

            processed_categories.append(category)

        HEADING.reset()
        HEADING.increase()
        text = build_rst(params)

        params_folder = Path(DOCS_SRC, "modules", category, "params")
        params_folder.mkdir(exist_ok=True, parents=True)
        with open(Path(params_folder, f"{module_name}.rst"), "w") as fout:
            fout.write(text)

        # relocate and title the module's apidoc stub
        process_module_file(category, module_name)

    # Generate general default parameters RST page
    HEADING.reset()
    HEADING.increase()
    general_defaults = Path(HADDOCK_REL_SRC_PATH, "modules", "defaults.yaml")
    general_params = read_from_yaml(general_defaults)
    text = build_rst(general_params)
    params_file = Path(DOCS_SRC, "modules", "general_module_params.rst")

    with open(params_file, "w") as fout:
        fout.write(text)

    # now libs, gear and core
    for folder in ("libs", "gear", "core"):
        # make directory if it does not exist
        target_path = Path(DOCS_SRC, "reference", folder)
        target_path.mkdir(exist_ok=True, parents=True)
        # collect rst files
        rst_files = DOCS_SRC.glob(f"haddock.{folder}*.rst")
        for rst_file in rst_files:
            target_rst = Path(target_path, rst_file.name)
            shutil.move(rst_file, target_rst)
            title_key = ".".join(rst_file.name.split(".")[1:-1])
            if title_key in REFERENCE_TITLE_DICT:
                title = REFERENCE_TITLE_DICT[title_key]
                change_title(target_rst, title)

    # Generate mandatory parameters RST page
    HEADING.reset()
    mandatory_defaults = Path(HADDOCK_REL_SRC_PATH, "core", "mandatory.yaml")
    mandatory_params = read_from_yaml(mandatory_defaults)

    for param in mandatory_params:
        mandatory_params[param]["default"] = (
            "No default assigned, this parameter is mandatory"
        )

    text = build_rst(mandatory_params)
    params_file = Path(DOCS_SRC, "reference", "core", "mandatory_parameters.rst")

    with open(params_file, "w") as fout:
        fout.write("Mandatory Parameters" + os.linesep)
        fout.write("====================" + os.linesep)
        fout.write(text)

    # now the command-line interfaces
    clients_folder = Path(DOCS_SRC, "clients")
    clients_folder.mkdir(exist_ok=True, parents=True)

    cli_rst_files = DOCS_SRC.glob("haddock.clis*.rst")
    for cli_rst_file in cli_rst_files:
        target_rst = Path(DOCS_SRC, "clients", cli_rst_file.name)
        shutil.move(cli_rst_file, target_rst)
        title_key = ".".join(cli_rst_file.name.split(".")[1:-1])
        if title_key in CLI_TITLE_DICT:
            title = CLI_TITLE_DICT[title_key]
            change_title(target_rst, title)


def do_text(name: str, param: ParamMap, level: str) -> str:
    """Create text from parameter dictionary."""
    text = [
        f"{name}",
        f"{level * len(name)}",
        "",
    ]
    text.append(f"| *default*: {param['default']!r}")
    text.append(f"| *type*: {param['type']}")
    text.append(f"| *title*: {param['title']}")

    # choices, min, max are not always present
    for key in ["choices", "min", "max"]:
        if key in param:
            text.append(f"| *{key}*: {param[key]}")

    text.append(f"| *short description*: {param['short']}")
    text.append(f"| *long description*: {param['long']}")
    text.append(f"| *group*: {param.get('group', 'No group assigned')}")
    text.append(f"| *explevel*: {param['explevel']}")
    text.append("")

    return os.linesep.join(text)


def loop_params(
    config: ParamMap, easy: list[str], expert: list[str], guru: list[str]
) -> tuple[list[str], list[str], list[str]]:
    """
    Treat parameters for module.

    *Important:* considers that some configuration files can have
    dictionaries with subparameters. However, there should NOT be more
    than one level of nesting in the configuration parameter files.
    """
    # sort parameters by name
    sorted_ = sorted(
        ((k, v) for k, v in config.items()),
        key=lambda x: x[0],
    )

    for name, data in sorted_:
        # case for nested parameters like `mol1` in topoaa
        if isinstance(data, Mapping) and "default" not in data:
            explevel = data["explevel"]
            new_title = [name, HEADING.current * len(name), ""]

            if explevel == "easy":
                easy.extend(new_title)
                sublist = easy
            elif explevel == "expert":
                expert.extend(new_title)
                sublist = expert
            elif explevel == "guru":
                guru.extend(new_title)
                sublist = guru
            elif explevel == "hidden":
                continue
            else:
                emsg = f"explevel {explevel!r} is not expected"
                raise AssertionError(emsg)

            data_text = (
                f"| *title*: {data['title']}",
                f"| *short description*: {data['short']}",
                f"| *long description*: {data['long']}",
                f"| *group*: {data['group']}",
                f"| *explevel*: {explevel}",
                "",
            )
            sublist.append(os.linesep.join(data_text))

            # create subparameters RST sorted by name
            data_sorted = sorted(
                ((k, v) for k, v in data.items()),
                key=lambda x: x[0],
            )
            for name2, param2 in data_sorted:
                if isinstance(param2, Mapping):
                    text = do_text(
                        f"{name}.{name2}",
                        param2,
                        level=HEADING.next,
                    )
                    sublist.append(text)

        # case for normal parameter
        elif isinstance(data, Mapping):
            explevel = data["explevel"]
            text = do_text(name, data, level=HEADING.current)

            if explevel == "easy":
                easy.append(text)
            elif explevel == "expert":
                expert.append(text)
            elif explevel == "guru":
                guru.append(text)
            elif explevel == "hidden":
                continue
            else:
                emsg = f"explevel {explevel!r} is not expected"
                raise AssertionError(emsg)
        else:
            emsg = f"Unexpected parameter behaviour: {name!r}"
            raise AssertionError(emsg)

    easy.append("")
    expert.append("")
    guru.append("")

    return easy, expert, guru


def build_rst(module_params: ParamMap) -> str:
    """Build .rst text."""
    easy = ["Easy", HEADING.current * 4, ""]
    expert = ["Expert", HEADING.current * 6, ""]
    guru = ["Guru", HEADING.current * 4, ""]

    HEADING.increase()
    easy, expert, guru = loop_params(module_params, easy, expert, guru)

    doc: list[str] = []
    for list_ in (easy, expert, guru):
        if len(list_) > 4:
            doc.extend(list_)

    text = os.linesep + os.linesep + os.linesep.join(doc)
    return text


if __name__ == "__main__":
    main()
