"""Greeting messages for the command line clients."""

import os
import random
import sys
from datetime import datetime
from functools import partial
from typing import Callable, Sequence

from haddock import contact_us, version


international_good_byes = [
    "Adéu-siau",
    "Agur",
    "Até logo",
    "Au revoir",
    "Ciao",
    "Dovidenia",
    "Good bye",
    "Tot ziens",
    "La revedere",
    "Adiós",
    "Do pobachennya",
    "Tchau",
    "再见",
]

# List of urls to be printed to the screen at the end of a workflow
# Do not hesitate to update / comment one of these
feedback_urls = {
    "GitHub issues": "https://github.com/haddocking/haddock3/issues",
    "BioExcel feedback": "https://www.bonvinlab.org/feedback",
    "BioExcel survey": "https://bioexcel.eu/bioexcel-survey-2025/",
    "BioExcel forum": "https://ask.bioexcel.eu/c/haddock/6",
}


DISCLAIMER = (
    "!! Some of the HADDOCK3 components use CNS (Crystallographic and NMR System)"
    f" which is free of use for non-profit applications. !!{os.linesep}"
    "!! For commercial use it is your own responsibility"
    f" to have a proper license. !!{os.linesep}"
    "!! For details refer to the DISCLAIMER file in the HADDOCK3 repository. !!"
)


def get_initial_greeting() -> str:
    """Create initial greeting message."""
    now = datetime.now().replace(second=0, microsecond=0)
    python_version = sys.version
    message = (
        f"""{os.linesep}"""
        f"""##############################################{os.linesep}"""
        f"""#                                            #{os.linesep}"""
        f"""#                 HADDOCK3                   #{os.linesep}"""
        f"""#                                            #{os.linesep}"""
        f"""##############################################{os.linesep}"""
        f"""{os.linesep}"""
        f"""{DISCLAIMER}{os.linesep}"""
        f"""{os.linesep}"""
        f"""Starting HADDOCK3 v{version} on {now}{os.linesep}"""
        f"""{os.linesep}"""
        f"""Python {python_version}{os.linesep}"""
    )
    return message


def get_greetings(
    options: Sequence[str], how_many: int = 3, sep: str = " ", exclamation: str = "!"
) -> str:
    """Get greeting messages."""
    n = how_many % len(options)
    return sep.join(s + exclamation for s in random.sample(options, k=n))


def get_adieu() -> str:
    """Create end-run greeting message."""
    end = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    bye = get_goodbye_greetings()
    message = f"Finished at {end}. {bye}"
    return message


def get_goodbye_help() -> str:
    """Create good-bye message with help."""
    end = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    bye = get_goodbye_greetings()
    message = f"Finished at {end}. For any help contact us at {contact_us}."
    message += f" {bye}."
    return message


def gen_feedback_messages(print_function: Callable) -> None:
    """Print list of feedbacks urls.

    Parameters
    ----------
    print_function : Callable
        The function used to print message on screen.
        This function must accept str as first argument.
    """
    print_function(
        (
            "Your feedback matters in Haddock3!"
            " Share your experience and help us grow:"
        )
    )
    for name, url in feedback_urls.items():
        print_function(f"{name}: {url}")


get_goodbye_greetings = partial(get_greetings, international_good_byes)
