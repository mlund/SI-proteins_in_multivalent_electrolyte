"""Execute the analysis notebooks in order, in place; never the MC simulation notebook."""

from pathlib import Path

import nbformat
from nbclient import NotebookClient

NOTEBOOKS = Path(__file__).resolve().parents[2] / "notebooks"
# 01-simulate compiles Faunus and runs for days; its output is already in data/processed.
ANALYSIS = ["02-pmf.ipynb"]


def analysis_notebooks():
    """Paths of the notebooks that reproduce figures from ``data/processed``."""
    return [NOTEBOOKS / name for name in ANALYSIS]


def main():
    """Run each analysis notebook from its own directory, saving outputs in place."""
    for path in analysis_notebooks():
        print(f"running {path.name}")
        nb = nbformat.read(path, as_version=4)
        NotebookClient(
            nb, timeout=600, resources={"metadata": {"path": str(NOTEBOOKS)}}
        ).execute()
        nbformat.write(nb, path)


if __name__ == "__main__":
    main()
