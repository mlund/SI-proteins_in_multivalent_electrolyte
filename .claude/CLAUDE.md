# Agent instructions

Supporting information for a scientific publication, archived on Zenodo.
Everything here must let a stranger reproduce the published results.

## Restructuring an existing project

When asked to turn an existing directory into supporting information:

1. Never rearrange the original in place. Clone this template into a new, separate directory
   (`git clone https://github.com/mlund/template-for-supporting-information <new-dir>`, then remove
   its `.git` and `git init`), and *copy* files into it. The original stays untouched as the reference.
2. Before copying, interview the user, one topic at a time, following the README checklist:
   - Title, description, keywords, and the associated publication (DOI, or journal if still in review).
   - Authors and contributors: full names, ORCIDs, affiliations, contributor roles.
   - Funding (`grants` in `.zenodo.json`).
   - License: code and data often warrant different licenses (e.g. MIT for code, CC-BY-4.0 for data,
     as `LICENSE` and `LICENSE-DATA`).
   - Structure: propose a mapping from the existing files onto the layout below and get it confirmed.
   - Documentation: what the reader needs to reproduce each result; data origin, formats and units.
   - Deposit size: see *Large files* below.
3. Copy, then fill in `.zenodo.json`, `CITATION.cff`, `pyproject.toml` and the README, and delete
   checklist items as they are done. Check that the notebooks run from the new directory.

## Large files

Decide per file type whether it belongs in git, and ask the user when unsure:

- Commit inputs, code, and small extracted features (e.g. per-frame observables, histograms, averages)
  that the figures are computed from. These are what most readers need.
- Keep large raw output (trajectories, checkpoints, binary dumps) out of git; add them to `.gitignore`.
  Ask whether they are needed at all, or whether extracted features plus the inputs to regenerate them
  suffice.

Zenodo deposits:

- Small repositories: use the GitHub-Zenodo integration; a GitHub release triggers a new version.
- Large data: upload archives manually, built by a script (e.g. `scripts/make_zenodo_zip.py`). A
  Zenodo record holds several files, and a new version can import unchanged files from the previous
  one, so split by how often things change: e.g. trajectories (never), extracted features (rarely),
  and code + results from one commit (often). Name each archive after a hash of its contents so an
  unchanged part keeps its name, include SHA256SUMS, and refuse to build from a dirty working tree.

## Skills

Before designing or writing code, load these skills (read each `SKILL.md` and follow it):

- TDD: https://github.com/mattpocock/skills/tree/main/skills/engineering/tdd
- Deep module design: https://github.com/mattpocock/skills/tree/main/skills/engineering/codebase-design

## Workflow

- Test-driven: write a failing test in `tests/`, make it pass, refactor. No logic without a test.
- Deep modules: small, simple public interfaces hiding substantial implementation.
  Prefer few deep modules over many shallow ones.
- Notebooks stay thin: they import from `src/analysis/`, call it and plot the results.
  Move any non-trivial logic out of notebooks and into tested code.

## Python

- Pythonic and idiomatic; target the `requires-python` in `pyproject.toml`.
- Prefer the standard library. Add a dependency only when it clearly pays for itself, and say why.
- Manage dependencies with uv only: `uv add <pkg>`, `uv add --dev <pkg>`. Never edit `uv.lock` by hand
  or use `pip`/conda. Commit `uv.lock`.
- OS-agnostic I/O: `pathlib.Path` for paths, explicit `encoding="utf-8"` on text I/O, paths relative
  to the repository, never absolute or user-specific.
- Public functions, classes and modules get numpy-style docstrings (ruff enforces this).
- Before finishing: `uv run ruff check --fix`, `uv run ruff format`, `uv run pytest`.

## Comments and docs

- Comments explain *why*, not *what*. If the code needs a *what* comment, rename or restructure it instead.
- Brief and terse, minimal jargon. Cite equations or papers where a choice follows from them.

## Layout

Keep directories ordered by the flow of the work; don't invent new top-level directories.

```
data/raw/        inputs, read-only, never modified by code
data/processed/  derived data, reproducible from raw + code
src/analysis/    importable, tested code
tests/           pytest tests mirroring src/analysis/
notebooks/       numbered in run order: 01-prepare.ipynb, 02-simulate.ipynb, ...
figures/         generated output, written by notebooks
```

## FAIR

- Use open, documented formats (CSV, JSON, HDF5, PDB...) over proprietary ones; state units.
- Record provenance: every file in `data/processed/` and `figures/` must be reproducible by a notebook.
- Keep `.zenodo.json` and `CITATION.cff` in sync when authors or metadata change.
