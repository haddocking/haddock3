"""CNS scripts util functions."""

import itertools
import math
from functools import partial
from os import linesep
from pathlib import Path

from haddock import EmptyPath, log
from haddock.core import cns_paths
from haddock.core.typing import Any, FilePath, FilePathT, Optional, Union
from haddock.libs import libpdb
from haddock.libs.libfunc import false, true
from haddock.libs.libmath import RandomNumberGenerator
from haddock.libs.libontology import PDBFile
from haddock.libs.libpdb import check_combination_chains
from haddock.libs.libutil import transform_to_list


RND = RandomNumberGenerator()


def generate_default_header(
    path: Optional[FilePath] = None,
) -> tuple[str, str, str, str, str, str]:
    """Generate CNS default header."""
    # TODO: Remove the `type: ignore` comments
    if path is not None:
        axis = load_axis(**cns_paths.get_axis(path))  # type: ignore
        link = load_link(Path(path, cns_paths.LINK_FILE))
        scatter = load_scatter(Path(path, cns_paths.SCATTER_LIB))
        tensor = load_tensor(**cns_paths.get_tensors(path))  # type: ignore
        trans_vec = load_trans_vectors(
            **cns_paths.get_translation_vectors(path)  # type: ignore
        )  # noqa: E501
        water_box = load_boxtyp20(cns_paths.get_water_box(path)["boxtyp20"])

    else:
        axis = load_axis(**cns_paths.axis)  # type: ignore
        link = load_link(cns_paths.link_file)
        scatter = load_scatter(cns_paths.scatter_lib)
        tensor = load_tensor(**cns_paths.tensors)  # type: ignore
        trans_vec = load_trans_vectors(**cns_paths.translation_vectors)  # type: ignore
        water_box = load_boxtyp20(cns_paths.water_box["boxtyp20"])

    return (
        link,
        trans_vec,
        tensor,
        scatter,
        axis,
        water_box,
    )


def find_desired_linkfiles(
        charged_nter: bool = False,
        charged_cter: bool = False,
        phosphate_5: bool = False,
        path: Optional[FilePath] = None,
        ) -> dict[str, Path]:
    """Find appropriate link files to use depending on terminis states.

    Parameters
    ----------
    charged_nter : bool, optional
        Must the Nter be charged ?, by default False
    charged_cter : bool, optional
        Must the Cter be charged ?, by default False
    phosphate_5 : bool, optional
        Must 5' be a phosphate ?, by default False
    path : Optional[FilePath], optional
        Path to where CNS topology/parameters are, by default None

    Returns
    -------
    linkfiles : dict[str, Path]
        Dict of CNS parameters/arguments/variable as keys
        and Path to link files to be used during topology
        generation.
    """
    # Set output variable
    linkfiles = {}
    # Logic to find appropriate link for proteins
    if charged_nter and charged_cter:
        prot_link_key = "NH3+,COO-"
    elif not charged_nter and charged_cter:
        prot_link_key = "NH,COO-"
    elif charged_nter and not charged_cter:
        prot_link_key= "NH3+,CO"
    elif not charged_nter and not charged_cter:
        prot_link_key = "NH,CO"
    # Point to corresponding file
    linkfiles["prot_link_infile"] = cns_paths.PROTEIN_LINK_FILES[prot_link_key]

    # Logic to find linkfile for dna
    nucl_link_key = "5'Phosphate" if phosphate_5 else "5'OH"
    # Point to corresponding file
    linkfiles["nucl_link_infile"] = cns_paths.NUCL_LINK_FILES[nucl_link_key]
    # Converts to real paths
    if path is not None:
        linkfiles = {key: Path(path, p) for key, p in linkfiles.items()}
    return linkfiles


def _is_nan(x: Any) -> bool:
    """Inspect if is nan."""
    try:
        return math.isnan(x)
    except (ValueError, TypeError):
        return False


def filter_empty_vars(v: Any) -> bool:
    """
    Filter empty variables.

    See: https://github.com/haddocking/haddock3/issues/162

    Returns
    -------
    bool
        Returns `True` if the variable is not empty, and `False` if
        the variable is empty. That is, `False` reflects those variables
        that should not be written in CNS.

    Raises
    ------
    TypeError
        If the type of `value` is not supported by CNS.
    """
    cases = (
        (lambda x: _is_nan(x), false),
        (lambda x: isinstance(x, str) and bool(x), true),
        (lambda x: isinstance(x, str) and not bool(x), false),
        (lambda x: isinstance(x, bool), true),  # it should return True
        (lambda x: isinstance(x, (EmptyPath, Path)), true),
        (lambda x: type(x) in (int, float), true),
        (lambda x: x is None, false),
    )

    for detect, give in cases:
        if detect(v):
            return give(v)
    else:
        emsg = f"Value {v!r} has a unknown type for CNS: {type(v)}."
        log.error(emsg)
        raise TypeError(emsg)


def load_workflow_params(
    param_header: str = f"{linesep}! Parameters{linesep}",
    **params: Any,
) -> str:
    """
    Write the values at the header section.

    "Empty variables" are ignored. These are defined accoring to
    :func:`filter_empty_vars`.

    Parameters
    ----------
    params : dict
        Dictionary containing the key:value pairs for the parameters to
        be written to CNS. Values cannot be of dictionary type.

    Returns
    -------
    param_header: str
        The string with the CNS parameters defined.
    """
    non_empty_parameters = (
        (k, v) for k, v in params.items()
        if filter_empty_vars(v)
        )

    # types besides the ones in the if-statements should not enter this loop
    for param, v in non_empty_parameters:
        param_header += write_eval_line(param, v)

    assert isinstance(param_header, str)
    return param_header


def write_eval_line(param: Any, value: Any, eval_line: str = "eval (${}={})") -> str:
    """Write the CNS eval line depending on the type of `value`."""
    eval_line += linesep

    if isinstance(value, bool):
        value = str(value).lower()
        return eval_line.format(param, value)

    elif isinstance(value, str):
        value = '"' + value + '"'
        return eval_line.format(param, value)

    elif isinstance(value, Path):
        value = '"' + str(value) + '"'
        return eval_line.format(param, value)

    elif isinstance(value, EmptyPath):
        return eval_line.format(param, '""')

    elif isinstance(value, (int, float)):
        return eval_line.format(param, value)

    else:
        emsg = f"Unexpected type when writing CNS header: {type(value)}"
        log.error(emsg)
        raise TypeError(emsg)


def load_link(mol_link: Path) -> str:
    """Add the link header."""
    return load_workflow_params(
        param_header=f"{linesep}! Link file{linesep}",
        prot_link_infile=mol_link,
    )


load_axis = partial(
    load_workflow_params, param_header=f"{linesep}! Axis{linesep}"
)  # noqa: E501
load_tensor = partial(
    load_workflow_params, param_header=f"{linesep}! Tensors{linesep}"
)  # noqa: E501
prepare_output = partial(
    load_workflow_params, param_header=f"{linesep}! Output structure{linesep}"
)  # noqa: E501
load_trans_vectors = partial(
    load_workflow_params, param_header=f"{linesep}! Translation vectors{linesep}"
)  # noqa: E501

load_ambig = partial(write_eval_line, "ambig_fname")
load_unambig = partial(write_eval_line, "unambig_fname")
load_hbond = partial(write_eval_line, "hbond_fname")
load_dihe = partial(write_eval_line, "dihe_f")
load_tensor_tbl = partial(write_eval_line, "tensor_tbl")


def load_scatter(scatter_lib: Path) -> str:
    """Add scatter library."""
    return load_workflow_params(
        param_header=f"{linesep}! Scatter lib{linesep}", scatter_lib=scatter_lib
    )


def load_boxtyp20(waterbox_param: Path) -> str:
    """Add boxtyp20 eval line."""
    return load_workflow_params(
        param_header=f"{linesep}! Water box{linesep}", boxtyp20=waterbox_param
    )


# This is used by docking
def prepare_multiple_input(
    pdb_input_list: list[str], 
    psf_input_list: list[str]
) -> str:
    """Prepare multiple input files."""
    input_str = f"{linesep}! Input structure{linesep}"
    for psf in psf_input_list:
        input_str += f"structure{linesep}"
        input_str += f"  @@{psf}{linesep}"
        input_str += f"end{linesep}"

    ncount = 1
    for pdb in pdb_input_list:
        input_str += f"coor @@{pdb}{linesep}"
        input_str += write_eval_line(f"input_pdb_filename_{ncount}", pdb)
        ncount += 1

    # check how many chains there are across all the PDBs
    chain_l: list[list[str]] = []
    for pdb in pdb_input_list:
        for element in libpdb.identify_chainseg(pdb):
            chain_l.append(element)
    ncomponents = len(set(itertools.chain(*chain_l)))
    input_str += write_eval_line("ncomponents", ncomponents)

    return input_str


# This is used by Topology and Scoring
def prepare_single_input(
    pdb_input: FilePath, psf_input: Union[None, FilePath, list[FilePathT]] = None
) -> str:
    """Input of the CNS file.

    This section will be written for any recipe even if some CNS variables
    are not used, it should not be an issue.
    """
    input_str = f"{linesep}! Input structure{linesep}"

    if psf_input:
        # if isinstance(psf_input, str):
        input_str += f"structure{linesep}"
        input_str += f"  @@{psf_input}{linesep}"
        input_str += f"end{linesep}"
        input_str += f"coor @@{pdb_input}{linesep}"
        if isinstance(psf_input, list):
            input_str += f"structure{linesep}"
            for psf in psf_input:
                input_str += f"  @@{psf}{linesep}"
            input_str += f"end{linesep}"

    # $file variable is still used by some CNS recipes, need refactoring!
    input_str += write_eval_line("file", pdb_input)
    segids, chains = libpdb.identify_chainseg(pdb_input)
    chainsegs = sorted(list(set(segids) | set(chains)))

    ncomponents = len(chainsegs)
    input_str += write_eval_line("ncomponents", ncomponents)

    for i, segid in enumerate(chainsegs, start=1):
        input_str += write_eval_line(f"prot_segid_{i}", segid)

    seed = RND.randint(100, 99999)
    input_str += write_eval_line("seed", seed)

    return input_str


def prepare_cns_input(
    model_number: int,
    input_element: Union[PDBFile, list[PDBFile]],
    step_path: FilePath,
    recipe_str: str,
    defaults: Any,
    identifier: str,
    ambig_fname: FilePath = "",
    native_segid: bool = False,
    cgtoaa: bool = False,
    default_params_path: Optional[Path] = None,
    debug: Optional[bool] = False,
    seed: Optional[int] = None,
) -> Union[Path, str]:
    """
    Generate the .inp file needed by the CNS engine.

    Parameters
    ----------
    model_number : int
        The number of the model. Will be used as file name suffix.

    input_element : `libs.libontology.Persisten`, list of those
    """
    # TODO: Refactor this function into smaller functions or classes
    # read the default parameters
    default_params = load_workflow_params(**defaults)
    default_params += write_eval_line("ambig_fname", ambig_fname)

    # write the PDBs
    pdb_list = [pdb.rel_path for pdb in transform_to_list(input_element)]

    # write the PSFs
    psf_list: list[Path] = []
    if isinstance(input_element, (list, tuple)):
        for pdb in input_element:
            if isinstance(pdb.topology, (list, tuple)):
                for psf in pdb.topology:
                    psf_fname = psf.rel_path
                    psf_list.append(psf_fname)
            else:
                if pdb.topology is None:
                    raise ValueError(f"Topology not found for pdb {pdb.rel_path}.")
                psf_fname = pdb.topology.rel_path
                psf_list.append(psf_fname)

    elif isinstance(input_element.topology, (list, tuple)):
        pdb = input_element  # for clarity
        if pdb.topology is None:
            raise ValueError(f"Topology not found for pdb {pdb.rel_path}.")
        for psf in pdb.topology:
            psf_fname = psf.rel_path
            psf_list.append(psf_fname)
    else:
        pdb = input_element  # for clarity
        if pdb.topology is None:
            raise ValueError(f"Topology not found for pdb {pdb.rel_path}.")
        psf_fname = pdb.topology.rel_path
        psf_list.append(psf_fname)

    aa_psf_list: list[Path] = []
    cgtoaa_tbl_list: list[Path] = []
    if cgtoaa==True:
        if isinstance(input_element.aa_topology, (list)):
            for psf in input_element.aa_topology:
                if psf is None:
                    raise ValueError(f"All-Atom Topology not found {input_element.rel_path}. "
                    "Conversion to all-atom requires a topology generated with [topoaa] and "
                    "[topocg].")
                else:
                    aa_psf_list.append(psf.rel_path.as_posix())
            for tbl in input_element.cgtoaa_tbl :
                if tbl is None:
                    raise ValueError(f"Coarse-Crain to All-Atom restraint file not found "
                    "{input_element.rel_path}. Conversion to all-atom requires a restraint file "
                    "generated with [topocg].")
                else:
                    cgtoaa_tbl_list.append(tbl.as_posix())
        else:
            pdb = input_element
            if pdb.aa_topology is None:
                raise ValueError(f"All-Atom Topology not found {input_element.rel_path}."
                "Conversion to all-atom requires a topology generated with [topoaa] and "
                "[topocg].")
            aa_psf_list.append(pdb.aa_topology.rel_path.as_posix())
            if pdb.cgtoaa_tbl is None:
                raise ValueError(f"Coarse-Crain to All-Atom restraint file not found for"
                " entry: {input_element.rel_path}. Conversion to all-atom requires a restraint file "
                "generated with [topocg].")
            cgtoaa_tbl_list.append(pdb.cgtoaa_tbl.as_posix())

    input_str = prepare_multiple_input(
        pdb_input_list=[str(p) for p in pdb_list],
        psf_input_list=[str(p) for p in psf_list],
    )
    
    if cgtoaa==True:
        for i in range(len(aa_psf_list)):
            # eval line for psf
            param_psf = "input_aa_psf_filename_" + str(i+1)
            input_str += write_eval_line(param_psf, aa_psf_list[i])
            # eval line for pdb
            param_pdb = "input_aa_pdb_filename_" + str(i+1)
            input_str += write_eval_line(param_pdb, aa_psf_list[i][:-4]+".pdb")
        for i in range(len(cgtoaa_tbl_list)):
            # eval line for tbl
            param_tbl = "input_cgtbl_filename_" + str(i+1)
            input_str += write_eval_line(param_tbl, cgtoaa_tbl_list[i])

    output_pdb_filename = f"{identifier}_{model_number}.pdb"

    output = f"{linesep}! Output structure{linesep}"
    output += write_eval_line("output_pdb_filename", output_pdb_filename)

    # prepare chain/seg IDs
    segid_str = ""
    if native_segid:
        if isinstance(input_element, (list, tuple)):
            chainid_list = check_combination_chains(input_element)

            for i, _chainseg in enumerate(chainid_list, start=1):
                segid_str += write_eval_line(f"prot_segid_{i}", _chainseg)

        else:
            chainid_list: list[str] = []
            segids, chains = libpdb.identify_chainseg(
                input_element.rel_path, sort=False
            )

            chainsegs = sorted(list(set(segids) | set(chains)))

            for i, _chainseg in enumerate(chainsegs, start=1):
                segid_str += write_eval_line(f"prot_segid_{i}", _chainseg)

    output += write_eval_line("count", model_number)

    if seed is None:
        seed = RND.randint(100, 99999)

    seed_str = write_eval_line("seed", seed)

    inp = default_params + input_str + seed_str + output + segid_str + recipe_str

    if not debug:
        return inp
    else:
        inp_file = Path(f"{identifier}_{model_number}.inp")
        inp_file.write_text(inp)
        return inp_file


def prepare_expected_pdb(
    model_obj: Union[PDBFile, tuple[PDBFile, ...]],
    model_nb: int,
    path: FilePath,
    identifier: str,
) -> PDBFile:
    """Prepare a PDBobject."""
    expected_pdb_fname = Path(path, f"{identifier}_{model_nb}.pdb")
    pdb = PDBFile(expected_pdb_fname, path=path)
    if isinstance(model_obj, tuple):
        pdb.topology = [p.topology for p in model_obj]
        pdb.aa_topology = [p.aa_topology for p in model_obj]
        pdb.cgtoaa_tbl = [p.cgtoaa_tbl for p in model_obj]
    else:
        pdb.topology = model_obj.topology
        pdb.seed = model_obj.seed
        pdb.aa_topology = model_obj.aa_topology
        pdb.cgtoaa_tbl = model_obj.cgtoaa_tbl
    return pdb
