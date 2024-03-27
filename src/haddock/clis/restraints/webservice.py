from base64 import b64decode
import io
import tempfile
from contextlib import redirect_stdout
from typing import Annotated, Union

from fastapi import FastAPI, Form, UploadFile
from pydantic import BaseModel, Field
from starlette.responses import PlainTextResponse
from uvicorn import run

from haddock.clis.restraints.passive_from_active import passive_from_active
from haddock.libs.librestraints import active_passive_to_ambig


"""
Webservice needs fastapi, python-multipart and uvicorn Python packages.
"""


def add_webservice_arguments(webservice_subcommand):
    webservice_subcommand.add_argument(
        "--host",
        help="Host for the webservice. Default is loopback.",
        default="127.0.0.1",
    )
    webservice_subcommand.add_argument(
        "-p",
        "--port",
        help="Port to run the webservice on. Default is 5000.",
        default=5000,
        type=int,
    )
    return webservice_subcommand


app = FastAPI()

# TODO add rate limit with slowapi package
# TODO if on Internet should have some authz


class PassiveFromActiveRequest(BaseModel):
    structure: str = Field(description="The structure file as base64 encoded string.")
    active: list[int] = Field(description="List of active restraints.")
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

        output = io.StringIO()
        with redirect_stdout(output):
            passive_from_active(
                structure=structure_file.name,
                active_list=",".join(map(str, request.active)),
                chain_id=request.chain,
                surface_list=",".join(map(str, request.surface)),
            )
        passive_as_string = output.getvalue()
        return passive_as_string.strip().split(" ")


class ActPassToAmbigRequest(BaseModel):
    active1: list[int] = Field(
        description="List of active residues for the first model."
    )
    active2: list[int] = Field(
        description="List of active residues for the second model."
    )
    passive1: list[int] = Field(
        description="List of passive residues for the first model."
    )
    passive2: list[int] = Field(
        description="List of passive residues for the second model."
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


def webservice(host, port):
    """Run the haddock3 restraints webservice."""
    run(app, port=port, host=host, log_level="debug")
