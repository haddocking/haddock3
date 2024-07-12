"""Webservice for the haddock3 restraints module.

Exposes the haddock3-restraints CLI subcommands as a RESTful webservice.
Also includes endpoint for PDB preprocessing.

Run with:

```shell
uvicorn --port 5000 haddock.clis.restraints.webservice:app
```

The Swagger UI is running at http://127.0.0.1:5000/docs .

To pass a PDB file in the JSON body of a request,
it needs to be gzipped and then base64 encoded.

A base64 encoded gzipped PDB file can made with:

```shell
cat examples/data/2oob.pdb | gzip | base64 -w 0 > 2oob.pdb.gz.base64
```

Background for PDB file handling:

To store a multiline string, like a pdb file,
in JSON we need to encode it in base64.
Base64 encoding make things 1.33 times bigger.
A pdb is text which can be compressed a lot.
To transfer less data we can compress the pdb
with gzip before base64 encoding.
For example the 2oob.pdb 74.8Kb becomes 101Kb when base64 encoded
while first gzip and then base64 encode it is 25.4Kb.
"""

import gzip
import io
import tempfile
from base64 import b64decode
from contextlib import redirect_stdout
from typing import Annotated

from fastapi import FastAPI, HTTPException, status
from fastapi.middleware.gzip import GZipMiddleware
from pdbtools.pdb_chain import alter_chain
from pdbtools.pdb_fixinsert import fix_insertions
from pdbtools.pdb_selaltloc import select_by_occupancy
from pdbtools.pdb_selchain import select_chain
from pdbtools.pdb_tidy import tidy_pdbfile
from pdbtools.pdb_delhetatm import remove_hetatm
from pdbtools.pdb_keepcoord import keep_coordinates
from pydantic import BaseModel, Field
from starlette.responses import PlainTextResponse

from haddock.clis.restraints.calc_accessibility import (
    apply_cutoff,
    get_accessibility,
)
from haddock.clis.restraints.restrain_bodies import (
    restrain_bodies as restrain_bodies_raw,
)
from haddock.libs.librestraints import (
    active_passive_to_ambig,
    check_parenthesis,
    passive_from_active_raw,
    validate_tbldata,
)


app = FastAPI()
app.add_middleware(GZipMiddleware, minimum_size=1000)


def unpacked_structure(
    structure: str,
) -> bytes:
    """Gunzips a base64 encoded string."""
    decoded = b64decode(structure)
    return gzip.decompress(decoded)


def unpacked_tbl(tbl: str) -> str:
    """Gunzips a base64 encoded tbl file contents."""
    return gzip.decompress(b64decode(tbl)).decode("utf-8")


Structure = Annotated[
    str,
    Field(
        description="The structure file as a base64 encoded gzipped string.",
        json_schema_extra=dict(
            contentMediaType="text/plain",
            contentEncoding="base64",
        )
    ),
]


class PassiveFromActiveRequest(BaseModel):
    structure: Structure
    active: list[int] = Field(
        description="List of active restraints.", examples=[[1, 2, 3]]
    )
    chain: str = Field(default="A", description="The chain identifier.")
    surface: list[int] = Field(default=[], description="List of surface restraints.")
    radius: float = Field(default=6.5, description="The radius from active.")


@app.post("/passive_from_active", tags=["restraints"])
def calculate_passive_from_active(
    request: PassiveFromActiveRequest,
) -> list[int]:
    """Calculate active restraints to passive restraints."""
    structure = unpacked_structure(request.structure)

    with tempfile.NamedTemporaryFile() as structure_file:
        structure_file.write(structure)

        try:
            passive = passive_from_active_raw(
                structure=structure_file.name,
                active=request.active,
                chain_id=request.chain,
                surface=request.surface,
                radius=request.radius,
            )
            return passive
        except Exception as e:
            raise HTTPException(
                status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e)
            ) from e


class ActPassToAmbigRequest(BaseModel):
    active1: list[int] = Field(
        description="List of active residues for the first model.", examples=[[1, 2, 3]]
    )
    active2: list[int] = Field(
        description="List of active residues for the second model.",
        examples=[[4, 5, 6]],
    )
    passive1: list[int] = Field(
        description="List of passive residues for the first model.",
        examples=[[7, 8, 9]],
    )
    passive2: list[int] = Field(
        description="List of passive residues for the second model.",
        examples=[[10, 11, 12]],
    )
    segid1: str = Field(default="A", description="Segid to use for the first model.")
    segid2: str = Field(default="B", description="Segid to use for the second model.")


@app.post("/actpass_to_ambig", response_class=PlainTextResponse, tags=["restraints"])
def calculate_actpass_to_ambig(
    request: ActPassToAmbigRequest,
) -> str:
    """Get the passive residues."""
    output = io.StringIO()
    with redirect_stdout(output):
        active_passive_to_ambig(
            active1=request.active1,
            passive1=request.passive1,
            active2=request.active2,
            passive2=request.passive2,
            segid1=request.segid1,
            segid2=request.segid2,
        )
    ambig = output.getvalue()
    return ambig.strip()


class RestrainBodiesRequest(BaseModel):
    structure: Structure
    exclude: list[str] = Field(
        default=[],
        description="Chains to exclude from the calculation.",
        examples=[["B"]],
    )


@app.post("/restrain_bodies", response_class=PlainTextResponse, tags=["restraints"])
def restrain_bodies(request: RestrainBodiesRequest) -> str:
    """Create distance restraints to lock several chains together."""
    structure = unpacked_structure(request.structure)
    with tempfile.NamedTemporaryFile() as structure_file:
        structure_file.write(structure)

        output = io.StringIO()
        with redirect_stdout(output):
            restrain_bodies_raw(
                structure=structure_file.name,
                exclude=request.exclude,
            )
        tbl = output.getvalue()
        return tbl.strip()


class CalcAccessibilityRequest(BaseModel):
    structure: Structure
    cutoff: float = Field(
        description="Relative cutoff for sidechain accessibility.", examples=[0.4]
    )


@app.post("/calc_accessibility", tags=["restraints"])
def calculate_accessibility(
    request: CalcAccessibilityRequest,
) -> dict[str, list[int]]:
    """Calculate the accessibility of the side chains and apply a cutoff."""
    structure = unpacked_structure(request.structure)

    with tempfile.NamedTemporaryFile() as structure_file:
        structure_file.write(structure)

        try:
            access_dic = get_accessibility(structure_file.name)
            # Filter residues based on accessibility cutoff
            result_dict = apply_cutoff(access_dic, request.cutoff)
            return result_dict
        except AssertionError as e:
            raise HTTPException(
                status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e)
            ) from e


class ValidateTblRequest(BaseModel):
    tbl: Annotated[
        str,
        Field(
            description="The TBL file as base64 encoded gzipped string.",
            json_schema_extra=dict(
                contentMediaType="text/plain",
                contentEncoding="base64",
            )
        ),
    ]
    pcs: bool = Field(
        default=False,
        description="Flag to indicate if the TBL file is in PCS mode.",
    )
    quick: bool = Field(
        default=False,
        description="Check global formatting before going line by line (opening/closing parenthesis and quotation marks.",
    )


@app.post("/validate_tbl", response_class=PlainTextResponse, tags=["restraints"])
def validate_tbl(
    request: ValidateTblRequest,
) -> str:
    tbl = unpacked_tbl(request.tbl)
    try:
        if request.quick:
            check_parenthesis(tbl)
        return validate_tbldata(tbl, request.pcs)
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e)
        ) from e


class PDBPreprocessRequest(BaseModel):
    structure: Structure
    from_chain: str = Field(description="Chains to keep", examples=["A"])
    to_chain: str = Field(description="New chain identifier", examples=["A"])
    delhetatm: bool = Field(description="Delete HETATM records", examples=[True], default=False)
    keepcoord: bool = Field(description="Remove all non-coordinate records", examples=[True], default=False)


@app.post("/preprocess_pdb", response_class=PlainTextResponse, tags=["pdb"])
def preprocess_pdb(request: PDBPreprocessRequest) -> str:
    """Preprocess a PDB file.

    Runs the following [pdbtools](http://www.bonvinlab.org/pdb-tools/) pipeline:

    ```shell
    cat pdb | pdb_tidy -strict | pdb_selchain -<from_chain> | pdb_chain -<to_chain> | pdb_fixinsert | pdb_selaltloc | pdb_tidy -strict
    ```

    or with `delhetatm` and `keepcoord` set to true:

    ```shell
    cat pdb | pdb_tidy -strict | pdb_selchain -<from_chain> | pdb_chain -<to_chain> | pdb_delhetatm | \
        pdb_fixinsert | pdb_keepcoord | pdb_selaltloc | pdb_tidy -strict
    ```

    """
    structure = unpacked_structure(request.structure).decode("latin_1")
    lines = structure.splitlines()
    lines = list(tidy_pdbfile(lines, strict=True))
    lines = list(select_chain(lines, request.from_chain))
    lines = list(alter_chain(lines, request.to_chain))
    if request.delhetatm:
        lines = list(remove_hetatm(lines))
    lines = list(fix_insertions(lines, []))
    if request.keepcoord:
        lines = list(keep_coordinates(lines))
    lines = list(select_by_occupancy(lines))
    lines = list(tidy_pdbfile(lines, strict=True))
    return "".join(lines)
