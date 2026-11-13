[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10437776.svg)](https://doi.org/10.5281/zenodo.10437776)
[![ci](https://github.com/mlund/SI-proteins_in_multivalent_electrolyte/actions/workflows/ci.yml/badge.svg)](https://github.com/mlund/SI-proteins_in_multivalent_electrolyte/actions/workflows/ci.yml)

# Supporting information: Anomalous protein–protein interactions in multivalent salt solution

C. Pasquier, M. Vazdar, J. Forsman, P. Jungwirth and M. Lund,
*J. Phys. Chem. B* **121**, 3000–3006 (2017),
[doi:10.1021/acs.jpcb.7b01051](https://doi.org/10.1021/acs.jpcb.7b01051).

Notebooks and tested code that reproduce the potentials of mean force (PMF) between two human
serum albumin (HSA) molecules in yttrium chloride solution, computed by Metropolis Monte Carlo
with [Faunus](https://github.com/mlund/faunus).

## Reproduce the figures

Install [uv](https://docs.astral.sh/uv/), then

``` bash
uv run analysis      # runs the analysis notebooks; writes figures/pmf-*.pdf
uv run jupyter lab   # or explore interactively
```

`uv` installs the Python version in `.python-version` and the exact packages in `uv.lock`.
`uv run analysis` does not rerun the Monte Carlo simulations, which take days; their output is in
`data/processed/pmf.csv`.

| Notebook | Produces | Run time |
|---|---|---|
| `notebooks/01-simulate.ipynb` | `data/processed/pmf.csv` (Faunus build, MC runs) | days; optional |
| `notebooks/02-pmf.ipynb` | `figures/pmf-0mM.pdf`, `pmf-100mM.pdf`, `pmf-200mM.pdf` | seconds |

## Layout

- `data/raw/faunus/` Monte Carlo inputs: `twobody.cpp` (Faunus program, 2017 API, commit
  `bdf6c7f7`) and `hsa-ph7.aam` (coarse-grained HSA, one site per amino acid at pH 7; columns: name,
  number, x, y, z (Å), charge (e), molar mass (g/mol), radius (Å)).
- `data/raw/gromacs/` GROMACS input files for the all-atom simulations of HSA.
- `data/processed/pmf.csv` PMFs, implicit 1:1 salt, specific Y³⁺–carboxylate attraction. Columns:
  `salt_M` added 1:1 salt (mol/l), `n_yttrium` number of Y³⁺ in the simulation cylinder (length
  563 Å, radius 50 Å; 1, 8, 16, 160 ions ≈ 0.4, 3, 6, 60 mM), `r_angstrom` protein mass-centre
  separation (Å), `pmf_kT` w(r)/k<sub>B</sub>T.
- `src/analysis/` Python code used by the notebooks: Faunus input (`faunus_input.py`), PMF
  analysis and plotting (`pmf.py`), notebook runner (`run.py`).
- `tests/` Tests, including the inputs and PMFs from v1.0.1 pinned as references.
- `notebooks/` Run in numbered order. `figures/` Generated figures.

## Development

``` bash
uv run pytest
uv run ruff check --fix && uv run ruff format
```

CI runs these and validates `.zenodo.json`.

## Changes since v1.0.1

The published numbers are unchanged; tests pin them against v1.0.1.

- The pandas pickle is now `data/processed/pmf.csv`; the notebook is split into simulation and
  analysis, with the logic in tested modules.
- Python environment with uv (replaces `environment.yml`); runs on current numpy, scipy and matplotlib.
- Adds the HSA structure `hsa-ph7.aam`, which v1.0.1 lacked.
- Corrects the article DOI and adds Zenodo and citation metadata.

## License

Code (`src/`, `tests/`, `notebooks/`, `data/raw/faunus/twobody.cpp`) is GPL-3.0-or-later ([LICENSE](LICENSE)).
Data, figures, simulation inputs and text are CC-BY-4.0 ([LICENSE-DATA](LICENSE-DATA)).
