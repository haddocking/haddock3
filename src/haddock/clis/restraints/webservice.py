"""Webservice for the haddock3 restraints module.

Webservice needs fastapi and uvicorn Python packages.

To run use:

uvicorn --port 5000 haddock.clis.restraints.webservice:app

Swagger UI at http://127.0.0.1:5000/docs

To base64 encode a file use:
base64 -w 0 file.txt
"""

from base64 import b64decode
import io
import tempfile
from contextlib import redirect_stdout

from fastapi import FastAPI, HTTPException, status
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

app = FastAPI(root_path="/restraints")

# TODO add rate limit with slowapi package
# TODO if on Internet should have some authz


class PassiveFromActiveRequest(BaseModel):
    structure: str = Field(
        description="The structure file as base64 encoded string.",
        contentMediaType="text/plain",
        contentEncoding="base64",
    )
    active: list[int] = Field(
        description="List of active restraints.", examples=[[1, 2, 3]]
    )
    chain: str = Field(default="A", description="The chain identifier.")
    surface: list[int] = Field(default=[], description="List of surface restraints.")


@app.post("/passive_from_active")
def calculate_passive_from_active(
    request: PassiveFromActiveRequest,
) -> list[int]:
    """
    Calculate active restraints to passive restraints.
    """
    with tempfile.NamedTemporaryFile() as structure_file:
        structure_file.write(b64decode(request.structure))

        try:
            passive = passive_from_active_raw(
                structure=structure_file.name,
                active=request.active,
                chain_id=request.chain,
                surface=request.surface,
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


@app.post("/actpass_to_ambig", response_class=PlainTextResponse)
def calculate_actpass_to_ambig(
    request: ActPassToAmbigRequest,
) -> str:
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
    structure: str = Field(
        description="The structure file as base64 encoded string.",
        contentMediaType="text/plain",
        contentEncoding="base64",
    )
    exclude: list[str] = Field(
        default=[],
        description="Chains to exclude from the calculation.",
        examples=[["B"]],
    )


@app.post("/restrain_bodies", response_class=PlainTextResponse)
def restrain_bodies(request: RestrainBodiesRequest) -> str:
    with tempfile.NamedTemporaryFile() as structure_file:
        structure_file.write(b64decode(request.structure))

        output = io.StringIO()
        with redirect_stdout(output):
            restrain_bodies_raw(
                structure=structure_file.name,
                exclude=request.exclude,
            )
        tbl = output.getvalue()
        return tbl.strip()


class CalcAccessibilityRequest(BaseModel):
    structure: str = Field(
        description="The structure file as base64 encoded string.",
        contentMediaType="text/plain",
        contentEncoding="base64",
    )
    cutoff: float = Field(
        description="Relative cutoff for sidechain accessibility.", examples=[0.4]
    )


@app.post("/calc_accessibility")
def calculate_accessibility(
    request: CalcAccessibilityRequest,
) -> dict[str, list[int]]:
    with tempfile.NamedTemporaryFile() as structure_file:
        structure_file.write(b64decode(request.structure))

        access_dic = get_accessibility(structure_file.name)
        # Filter residues based on accessibility cutoff
        result_dict = apply_cutoff(access_dic, request.cutoff)
        return result_dict


class ValidateTblRequest(BaseModel):
    tbl: str = Field(
        description="The TBL file as base64 encoded string.",
        contentMediaType="text/plain",
        contentEncoding="base64",
    )
    pcs: bool = Field(
        default=False,
        description="Flag to indicate if the TBL file is in PCS mode.",
    )
    quick: bool = Field(
        default=False,
        description="Check global formatting before going line by line (opening/closing parenthesis and quotation marks.",
    )


@app.post("/validate_tbl", response_class=PlainTextResponse)
def validate_tbl(
    request: ValidateTblRequest,
) -> str:
    tbl = b64decode(request.tbl).decode("utf-8")
    try:
        if request.quick:
            check_parenthesis(tbl)
        return validate_tbldata(tbl, request.pcs)
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e)
        ) from e
