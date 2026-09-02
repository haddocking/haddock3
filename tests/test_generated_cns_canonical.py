"""Golden canonical forms built from real HADDOCK CNS inputs."""

import os
from collections import Counter
from copy import deepcopy
from difflib import SequenceMatcher
from pathlib import Path

import pytest

import haddock
from haddock.libs.libcns import (
    derive_seed,
    find_desired_linkfiles,
    prepare_cns_input,
)
from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.libs.libseamless import canonical_mapping_for_job
from haddock.libs.libsubprocess import CNSJob
from haddock.modules.refinement.cgtoaa import HaddockModule as Cgtoaa
from haddock.modules.refinement.emref import HaddockModule as Emref
from haddock.modules.refinement.flexref import HaddockModule as Flexref
from haddock.modules.refinement.mdref import HaddockModule as Mdref
from haddock.modules.sampling.rigidbody import HaddockModule as Rigidbody
from haddock.modules.scoring.emscoring import HaddockModule as Emscoring
from haddock.modules.scoring.mdscoring import HaddockModule as Mdscoring
from haddock.modules.topology.topoaa import HaddockModule as Topoaa
from haddock.modules.topology.topoaa import generate_topology as generate_topoaa
from haddock.modules.topology.topocg import HaddockModule as Topocg
from haddock.modules.topology.topocg import generate_topology as generate_topocg

from . import golden_data


SOURCE_GOLDEN_DATA = Path(golden_data).resolve()
GOLDEN_DIR = SOURCE_GOLDEN_DATA / "cns_canonical"
#: The shapes whose recipes read ``$seed``.
SEEDED_SHAPES = frozenset(
    {"rigidbody", "flexref", "emref", "mdref", "mdscoring"}
)
GENERIC_MODULES = {
    "rigidbody": Rigidbody,
    "flexref": Flexref,
    "emref": Emref,
    "mdref": Mdref,
    "emscoring": Emscoring,
    "mdscoring": Mdscoring,
}


@pytest.mark.parametrize(
    ("module_class", "constructed_parameters"),
    (
        (
            Rigidbody,
            {
                "int_1_2",
                "nrair_1",
                "rair_sta_1_1",
                "rair_end_1_1",
                "c6sym_seg6_1",
            },
        ),
        (
            Flexref,
            {
                "int_1_2",
                "seg_sta_1_1",
                "seg_end_1_1",
                "nseg1",
                "ncs_sta1_1",
                "c6sym_seg6_1",
            },
        ),
        (
            Emref,
            {
                "int_1_2",
                "seg_sta_1_1",
                "seg_end_1_1",
                "nseg1",
                "ncs_sta1_1",
                "c6sym_seg6_1",
            },
        ),
        (
            Mdref,
            {
                "int_1_2",
                "seg_sta_1_1",
                "seg_end_1_1",
                "nseg1",
                "ncs_sta1_1",
                "c6sym_seg6_1",
            },
        ),
        (
            Cgtoaa,
            {"int_1_2", "seg_sta_1_1", "seg_end_1_1", "nseg1"},
        ),
    ),
)
def test_constructed_cns_parameter_families_are_included(
    module_class,
    constructed_parameters,
    tmp_path,
):
    """Every affected module binds representative spliced CNS symbols."""
    module = module_class(0, tmp_path)

    assert constructed_parameters <= module.cns_params().keys()


@pytest.mark.parametrize(
    "shape",
    (
        "topoaa",
        "topocg",
        "rigidbody",
        "flexref",
        "emref",
        "mdref",
        "emscoring",
        "mdscoring",
        "cgtoaa",
    ),
)
def test_generated_cns_input_matches_canonical_golden(shape, tmp_path, monkeypatch):
    """Canonicalize a generated production input and compare its canonical form."""
    input_path = _stage_inputs(tmp_path)
    work_path = tmp_path / "step"
    work_path.mkdir()
    monkeypatch.chdir(work_path)
    mapping, script, module = _generated_mapping(shape, work_path, input_path)
    observed = _canonical_golden_text(mapping, script, module)
    golden = GOLDEN_DIR / f"{shape}.canonical"
    if os.environ.get("HADDOCK_UPDATE_CNS_GOLDENS") == "1":
        golden.parent.mkdir(parents=True, exist_ok=True)
        golden.write_text(observed, encoding="utf-8")

    assert observed == golden.read_text(encoding="utf-8")


#: What a golden form is, for whoever reads one in a review.
GOLDEN_PREAMBLE = """\
# The canonical form of one real generated CNS job, in reviewable pieces.
#
# [pins] binds each canonical name to the file that occupies it, by basename
# and content checksum.  Both halves are identity: which pin an input occupies
# and what that input contains.  A change to a file a recipe includes moves a
# checksum here, which is the only place it would otherwise be invisible.
#
# [outputs] is the declared output shape.
#
# [recipe rewrites] is everything canonicalization does to the module's CNS
# recipe, as before/after pairs with the number of places each occurs.  The
# recipe itself is not copied here -- it is in the tree already, and copying it
# would mirror every edit to it into a file that has no opinion about the
# change.  Editing a recipe therefore leaves this section alone unless it
# changes what canonicalization has to do.
#
# [canonical parameter header] is the generated part of the input -- the
# parameter set, the seed, the pin references, the erasures -- which exists
# nowhere else, so it is kept verbatim.
#
# Regenerate with HADDOCK_UPDATE_CNS_GOLDENS=1; the diff is the review.
"""


def _canonical_golden_text(mapping, script: str, module) -> str:
    """Render the canonical form of one generated job as a reviewable file."""
    canonical_lines = mapping.canonical_script.splitlines(keepends=True)
    generated_lines = script.splitlines(keepends=True)
    recipe_lines = module.recipe_str.splitlines(keepends=True)

    # A module recipe is spliced into the generated input unchanged, as its
    # tail.  That is what makes it separable from the generated part at all,
    # so it is asserted rather than assumed.
    header_length = len(generated_lines) - len(recipe_lines)
    assert generated_lines[header_length:] == recipe_lines

    boundary, rewrites = _recipe_rewrites(
        generated_lines, canonical_lines, header_length
    )

    sections = [GOLDEN_PREAMBLE, "\n[pins]\n"]
    for dependency in sorted(
        mapping.dependencies, key=lambda dependency: dependency.canonical_name
    ):
        # An install-tree pin is named after the file already; only a job
        # input needs to say which file it stands for.
        name = dependency.original_path.name
        binding = dependency.canonical_name
        if not binding.endswith(name):
            binding = f"{binding} <- {name}"
        sections.append(f"{binding}  {dependency.checksum}\n")
    sections.append(f"\n[outputs] {mapping.output_shape}\n")
    sections.extend(f"{name}\n" for name in mapping.canonical_output_names)
    sections.append(f"\n[recipe rewrites] {_recipe_label(module)}\n")
    if not rewrites:
        sections.append("(none)\n")
    for (before, after), count in sorted(rewrites.items()):
        source, canonical = before.rstrip("\n"), after.rstrip("\n")
        sections.append(f"{count}x  - {source}\n")
        sections.append(f"    + {canonical}\n")
    sections.append("\n[canonical parameter header]\n")
    sections.append("".join(canonical_lines[:boundary]))
    return "".join(sections)


def _recipe_rewrites(generated_lines, canonical_lines, header_length):
    """Find the recipe in the canonical script and count what changed in it.

    Canonicalization drops whole assignments out of the generated header, so
    the two scripts do not line up by index and are aligned instead.  Within
    the recipe it only rewrites lines in place; a rule that started adding or
    removing lines there would misalign these sections silently, so that is
    asserted too.
    """
    matcher = SequenceMatcher(None, generated_lines, canonical_lines, autojunk=False)
    boundary = None
    rewrites: Counter = Counter()
    for tag, source, source_end, canonical, canonical_end in matcher.get_opcodes():
        if source_end <= header_length:
            continue
        if source < header_length:
            assert tag == "equal", "the recipe boundary is inside a rewritten block"
            boundary = canonical + (header_length - source)
            continue
        if boundary is None:
            boundary = canonical
        if tag == "equal":
            continue
        assert tag == "replace" and source_end - source == canonical_end - canonical, (
            "canonicalization added or removed a line inside the recipe"
        )
        rewrites.update(
            zip(
                generated_lines[source:source_end],
                canonical_lines[canonical:canonical_end],
            )
        )
    assert boundary is not None
    return boundary, rewrites


def _recipe_label(module) -> str:
    """Name a module recipe by its place in the package, not on this machine."""
    package_root = Path(haddock.__file__).resolve().parent.parent
    return str(Path(module.cns_protocol_path).resolve().relative_to(package_root))


def _stage_inputs(tmp_path: Path) -> Path:
    """Create the sibling input directory used by production step paths."""
    input_path = tmp_path / "inputs"
    input_path.mkdir()
    sources = {
        name: SOURCE_GOLDEN_DATA / name
        for name in (
            "e2aP_1F3G_haddock.pdb",
            "e2aP_1F3G_haddock.psf",
            "hpr_ensemble_1_haddock.pdb",
            "hpr_ensemble_1_haddock.psf",
            "example_ambig_1.tbl",
        )
    }
    integration_data = Path(__file__).parents[1] / "integration_tests" / "golden_data"
    sources.update(
        {
            name: integration_data / name
            for name in ("e2a_haddock_cg.pdb", "e2a_haddock_cg.psf")
        }
    )
    for name, source in sources.items():
        (input_path / name).symlink_to(source.resolve())
    return input_path


def _generated_mapping(shape: str, work_path: Path, input_path: Path):
    if shape == "topoaa":
        module = Topoaa(0, work_path)
        input_pdb = input_path / "e2aP_1F3G_haddock.pdb"
        mol_params = deepcopy(module.params["mol1"])
        charged_nter = mol_params.pop("charged_nter")
        charged_cter = mol_params.pop("charged_cter")
        phosphate_5 = mol_params.pop("5_phosphate")
        mol_params.update(
            find_desired_linkfiles(
                charged_nter,
                charged_cter,
                phosphate_5,
                module.toppar_path,
            )
        )
        script = generate_topoaa(
            input_pdb,
            module.recipe_str,
            module.cns_params(),
            mol_params,
            default_params_path=module.toppar_path,
            write_to_disk=False,
        )
        output_stem = f"{input_pdb.stem}_haddock"
        outputs = [Path(f"{output_stem}.pdb"), Path(f"{output_stem}.psf")]
    elif shape == "topocg":
        module = Topocg(0, work_path)
        input_pdb = input_path / "e2aP_1F3G_haddock.pdb"
        script = generate_topocg(
            input_pdb,
            str(work_path),
            module.recipe_str,
            module.cns_params(),
            module.params["mol1"],
            default_params_path=module.toppar_path,
            write_to_disk=False,
            shape=True,
        )
        outputs = [Path(input_pdb.name), Path(f"{input_pdb.stem}.psf")]
    else:
        module, input_element, cgtoaa = _generic_input(
            shape,
            work_path,
            input_path,
        )
        # The seeded shapes are given the seed production derives for them,
        # so the golden form pins the derivation as well as the layout.
        seed = (
            derive_seed(module.params["iniseed"], input_element)
            if shape in SEEDED_SHAPES
            else None
        )
        script = prepare_cns_input(
            1,
            input_element,
            work_path,
            module.recipe_str,
            module.cns_params(),
            shape,
            default_params_path=module.toppar_path,
            native_segid=shape == "rigidbody",
            cgtoaa=cgtoaa,
            seed=seed,
        )
        outputs = [Path(f"{shape}_1.pdb")]

    job = CNSJob(
        script,
        envvars=module.default_envvars(),
        output_files=outputs,
    )
    return canonical_mapping_for_job(job), script, module


def _generic_input(shape: str, work_path: Path, input_path: Path):
    if shape == "cgtoaa":
        module = Cgtoaa(0, work_path)
        model = PDBFile(
            file_name="e2a_haddock_cg.pdb",
            path=input_path,
            topology=Persistent(
                file_name="e2a_haddock_cg.psf",
                path=input_path,
                file_type=Format.TOPOLOGY,
            ),
        )
        model.aa_topology = Persistent(
            file_name="e2aP_1F3G_haddock.psf",
            path=input_path,
            file_type=Format.TOPOLOGY,
        )
        model.cgtoaa_tbl = (input_path / "example_ambig_1.tbl").resolve()
        return module, model, True

    module = GENERIC_MODULES[shape](0, work_path)
    inputs = [
        _model("e2aP_1F3G_haddock", input_path),
        _model("hpr_ensemble_1_haddock", input_path),
    ]
    return module, inputs, False


def _model(stem: str, input_path: Path) -> PDBFile:
    return PDBFile(
        file_name=f"{stem}.pdb",
        path=input_path,
        topology=Persistent(
            file_name=f"{stem}.psf",
            path=input_path,
            file_type=Format.TOPOLOGY,
        ),
    )
