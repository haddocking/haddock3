"""
Run paired rigidbody vs gdock examples and compare their final caprieval
results.

USAGE:

    $ python devtools/compare_gdock_rigidbody.py --list
    $ python devtools/compare_gdock_rigidbody.py [-p PAIR ...] [--compare-only] [--ncores N]

With no `-p` given, all pairs are processed.

Both configs in a pair are run with `mode = "local"` and `ncores = N`
(default 8), regardless of what is set in the original config files, so
that the comparison is run under identical conditions.

Requires the haddock3 python environment to be active.
"""

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


EXAMPLES_DIR = Path(__file__).resolve().parent.parent / "examples"


@dataclass(frozen=True)
class Pair:
    """A rigidbody/gdock pair of example configs to compare."""

    name: str
    folder: str
    rigid_cfg: str
    gdock_cfg: str


PAIRS = [
    Pair(
        name="CDR-accessible",
        folder="docking-antibody-antigen",
        rigid_cfg="docking-antibody-antigen-CDR-accessible-full.cfg",
        gdock_cfg="CDR-accessible-gdock.toml",
    ),
    Pair(
        name="CDR-NMR-CSP",
        folder="docking-antibody-antigen",
        rigid_cfg="docking-antibody-antigen-CDR-NMR-CSP-full.cfg",
        gdock_cfg="CDR-NMR-CSP-gdock.toml",
    ),
    Pair(
        name="protein-protein",
        folder="docking-protein-protein",
        rigid_cfg="docking-protein-protein-full.cfg",
        gdock_cfg="protein-protein-gdock.toml",
    ),
    Pair(
        name="protein-peptide",
        folder="docking-protein-peptide",
        rigid_cfg="docking-protein-peptide-full.cfg",
        gdock_cfg="protein-peptide-gdock.toml",
    ),
    Pair(
        name="Para-Epi",
        folder="docking-nanobody-antigen",
        rigid_cfg="docking-nanobody-antigen-Para-Epi-full.cfg",
        gdock_cfg="Para-Epi-gdock.toml",
    ),
    Pair(
        name="CDR-loose-epi",
        folder="docking-nanobody-antigen",
        rigid_cfg="docking-nanobody-antigen-CDR-loose-epi-full.cfg",
        gdock_cfg="CDR-loose-epi-gdock.toml",
    ),
    Pair(
        name="CDR-mutagenesis-epi",
        folder="docking-nanobody-antigen",
        rigid_cfg="docking-nanobody-antigen-CDR-mutagenesis-epi-full.cfg",
        gdock_cfg="CDR-mutagenesis-epi-gdock.toml",
    ),
]


def run_dir_of(cfg_path):
    """Extract the `run_dir = "..."` value from a config file."""
    text = cfg_path.read_text()
    match = re.search(r'^\s*run_dir\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if not match:
        raise ValueError(f"Could not find run_dir in {cfg_path}")
    return match.group(1)


def set_mode_and_ncores(text, ncores):
    """Force `mode = "local"` and `ncores = <ncores>` in a config's text."""
    if re.search(r"^\s*mode\s*=", text, re.MULTILINE):
        text = re.sub(
            r"^\s*mode\s*=.*$", 'mode = "local"', text, count=1, flags=re.MULTILINE
        )
    else:
        text = f'mode = "local"\n{text}'

    if re.search(r"^\s*ncores\s*=", text, re.MULTILINE):
        text = re.sub(
            r"^\s*ncores\s*=.*$",
            f"ncores = {ncores}",
            text,
            count=1,
            flags=re.MULTILINE,
        )
    else:
        text = f"ncores = {ncores}\n{text}"

    return text


def run_cfg(folder, cfg, ncores):
    """Remove any previous run directory and execute `haddock3 <cfg>`.

    A temporary copy of the config is used with `mode = "local"` and
    `ncores = <ncores>` enforced, so that rigidbody and gdock runs are
    executed under identical conditions.
    """
    cfg_path = folder / cfg
    run_dir = folder / run_dir_of(cfg_path)

    text = set_mode_and_ncores(cfg_path.read_text(), ncores)
    tmp_cfg = folder / f".tmp_{cfg}"
    tmp_cfg.write_text(text)

    print(
        f">>> [{folder.name}] running {cfg} (run_dir={run_dir.name}, ncores={ncores})"
    )
    shutil.rmtree(run_dir, ignore_errors=True)
    try:
        subprocess.run(["haddock3", tmp_cfg.name], cwd=folder, check=True)
    finally:
        tmp_cfg.unlink()


def last_caprieval_dir(run_dir):
    """Return the highest-numbered `*_caprieval` directory in run_dir."""
    candidates = sorted(run_dir.glob("*_caprieval"))
    return candidates[-1] if candidates else None


def print_top_cluster(folder, cfg, label):
    """Print the cluster_rank=1 row of capri_clt.tsv from the last caprieval step."""
    run_dir = folder / run_dir_of(folder / cfg)
    capri_dir = last_caprieval_dir(run_dir)

    if capri_dir is None:
        print(f"{label}: no caprieval output found in {run_dir}")
        return

    print(f"--- {label} ({capri_dir.relative_to(folder)}) ---")
    clt_tsv = capri_dir / "capri_clt.tsv"
    if not clt_tsv.exists():
        print("no capri_clt.tsv found")
        return

    lines = clt_tsv.read_text().splitlines()
    header_idx = next(
        (i for i, line in enumerate(lines) if line.startswith("cluster_rank")), None
    )
    if header_idx is None or header_idx + 1 >= len(lines):
        print("no cluster rows found in capri_clt.tsv")
        return

    header = lines[header_idx].split("\t")
    values = lines[header_idx + 1].split("\t")
    for key, value in zip(header, values):
        print(f"  {key}: {value}")


def compare_pair(folder, rigid_cfg, gdock_cfg, name):
    print()
    print(f"==== {name} ({folder.name}) ====")
    print_top_cluster(folder, rigid_cfg, "rigidbody")
    print_top_cluster(folder, gdock_cfg, "gdock")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument(
        "-p",
        "--pair",
        action="append",
        dest="pairs",
        metavar="NAME",
        help="Name of a pair to process (can be given multiple times). "
        "If omitted, all pairs are processed.",
    )
    ap.add_argument(
        "-l",
        "--list",
        action="store_true",
        help="List the available pairs and exit.",
    )
    ap.add_argument(
        "--compare-only",
        action="store_true",
        help="Only compare existing run directories, do not run anything.",
    )
    ap.add_argument(
        "--ncores",
        type=int,
        default=8,
        help="Number of cores to use for both runs (default: 8). "
        'Also forces mode = "local" on both configs.',
    )
    args = ap.parse_args()

    if args.list:
        for pair in PAIRS:
            print(f"  - {pair.name}")
        return

    selected = set(args.pairs) if args.pairs else None

    for pair in PAIRS:
        if selected is not None and pair.name not in selected:
            continue

        folder = EXAMPLES_DIR / pair.folder

        if not args.compare_only:
            run_cfg(folder, pair.rigid_cfg, args.ncores)
            run_cfg(folder, pair.gdock_cfg, args.ncores)

        compare_pair(folder, pair.rigid_cfg, pair.gdock_cfg, pair.name)


if __name__ == "__main__":
    sys.exit(main())
