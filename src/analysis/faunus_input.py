"""Faunus ``hsa.json`` input for two HSA molecules with Y³⁺ in a periodic cylinder.

Constants and arithmetic mirror the 2017 notebook exactly, so regenerated inputs are
bit-identical to the ones behind the published PMFs (pinned in ``tests/reference``).
"""

import math

PH = 7.0
TEMPERATURE = 298.15  # K
KT = (1.3806e-23) * (6.022e23) * TEMPERATURE / 1000  # kJ/mol
AVOGADRO = 6.02214129e23  # 2017 CODATA value used for the published runs
CYLINDER_LENGTH = 563  # Å
CYLINDER_RADIUS = 50  # Å
VOLUME_L = (math.pi * (CYLINDER_RADIUS**2) * CYLINDER_LENGTH) * 1e-27
COUNTERIONS = 16  # Na⁺ neutralizing two HSA at pH 7 (net charge −8 each)

EPS = 0.05 * KT  # LJ well depth for all sites (kJ/mol)
EPS_ION = EPS / 10
EPS_Y_CARBOXYLATE = 1.6 * KT  # specific Y³⁺–carboxylate attraction

# name: (charge, radius Å, molar mass g/mol); None = no mass given in 2017 input
_SITES = {
    "H3PO4": (0, 2.0, None),
    "H2PO4": (-1, 2.0, None),
    "HPO4": (-2, 2.0, None),
    "PO4": (-3, 2.0, None),
    "BPTI": (7.3, 12.29, None),
    "I": (-1, 2.0, 1),
    "SCN": (-1, 2.0, 1),
    "ASP": (-1, 3.6, 110),
    "HASP": (0, 3.6, 110),
    "LASP": (2, 3.6, 110),
    "CTR": (-1, 2.0, 16),
    "HCTR": (0, 2.0, 16),
    "GLU": (-1, 3.8, 122),
    "HGLU": (0, 3.8, 122),
    "LGLU": (2, 3.8, 122),
    "HIS": (0, 3.9, 130),
    "HHIS": (1, 3.9, 130),
    "NTR": (0, 2.0, 14),
    "HNTR": (1, 2.0, 14),
    "TYR": (-1, 4.1, 154),
    "HTYR": (0, 4.1, 154),
    "LYS": (0, 3.7, 116),
    "HLYS": (1, 3.7, 116),
    "CYS": (0, 3.6, 103),
    "CYF": (-1, 3.6, 103),
    "HCYS": (0, 3.6, 103),
    "ARG": (0, 4.0, 144),
    "HARG": (1, 4.0, 144),
    "ALA": (0, 3.1, 66),
    "ILE": (0, 3.6, 102),
    "LEU": (0, 3.6, 102),
    "MET": (0, 3.8, 122),
    "PHE": (0, 3.9, 138),
    "PRO": (0, 3.4, 90),
    "TRP": (0, 4.3, 176),
    "VAL": (0, 3.4, 90),
    "SER": (0, 3.3, 82),
    "THR": (0, 3.5, 94),
    "ASN": (0, 3.6, 108),
    "GLN": (0, 3.8, 120),
    "GLY": (0, 2.9, 54),
    "SPH": (1, 3.0, 30),
    "CA": (2, 2.12, 40.08),
    "CYT": (0, 3.6, 103),
}

# (bound, free, pKa)
_TITRATIONS = {
    "H-Asp": ("HASP", "ASP", 4.0),
    "H-Ctr": ("HCTR", "CTR", 2.6),
    "H-Glu": ("HGLU", "GLU", 4.4),
    "H-His": ("HHIS", "HIS", 6.3),
    "H-Arg": ("HARG", "ARG", 12.0),
    "H-Ntr": ("HNTR", "NTR", 7.5),
    "H-Cys": ("HCYS", "CYS", 10.8),
    "H-Tyr": ("HTYR", "TYR", 9.6),
    "H-Lys": ("HLYS", "LYS", 10.4),
    "K1": ("H3PO4", "H2PO4", 2.12),
    "K2": ("H2PO4", "HPO4", 7.21),
    "K3": ("HPO4", "PO4", 12.67),
}


def _displacements(n_yttrium, explicit_ions):
    """Return (protein dp, protein rotation, Na dp, Y/Cl dp), tuned for acceptance."""
    if explicit_ions:
        return (80, 5, 80, 50) if 1 <= n_yttrium < 10 else (20, 3, 80, 30)
    if 1 <= n_yttrium < 35:
        return 80, 6, 70, 70
    if 35 <= n_yttrium < 60:
        return 40, 3, 70, 70
    return 20, 1, 70, 70


def _ionic_strength(n_yttrium, salt_M, explicit_ions):
    """Debye–Hückel ionic strength (mol/l) from the ions not simulated explicitly."""
    if explicit_ions:
        return 1e-20  # all ions explicit; Faunus needs a finite value
    # Monovalent Cl⁻ from YCl₃ plus Na⁺ counterions; Y³⁺ itself is explicit.
    monovalent = ((3 * n_yttrium + COUNTERIONS) / AVOGADRO) / VOLUME_L
    ionic = 0.5 * monovalent
    ionic += salt_M
    return ionic


def _atomlist(dp_na, dp_yt):
    atoms = {}
    for name, (q, r, mw) in _SITES.items():
        atoms[name] = {"eps": EPS, "q": q, "r": r} | (
            {"mw": mw} if mw is not None else {}
        )
    atoms["Na"] = {"eps": EPS_ION, "q": 1, "r": 1.9, "mw": 22.99, "dp": dp_na}
    atoms["Cl"] = {"eps": EPS_ION, "q": -1, "r": 1.7, "mw": 35.45, "dp": dp_yt}
    atoms["Yt"] = {"eps": EPS_ION, "q": 3, "r": 0.0, "mw": 88, "dp": dp_yt}
    return atoms


def faunus_input(
    n_yttrium,
    salt_M,
    *,
    explicit_ions,
    specific,
    micro,
    structure="hsa-ph7.aam",
):
    """Build the Faunus ``twobody`` input for one condition.

    Parameters
    ----------
    n_yttrium : int
        Number of Y³⁺ ions in the cylinder.
    salt_M : float
        Added 1:1 salt (mol/l), treated implicitly via the Debye length.
        Ignored when `explicit_ions` is true, as in the published runs.
    explicit_ions : bool
        Simulate Cl⁻ (from YCl₃) and Na⁺ counterions explicitly instead of
        through the ionic strength.
    specific : bool
        Add the specific Y³⁺–carboxylate (ASP, GLU, C-terminus) attraction.
    micro : int
        Monte Carlo steps per macro step (10 macro steps).
    structure : str
        Path to the coarse-grained HSA structure, relative to the run directory.

    Returns
    -------
    dict
        JSON-serializable input for ``hsa.json``.
    """
    dpp, drot, dp_na, dp_yt = _displacements(n_yttrium, explicit_ions)
    protein = {"structure": structure, "Ninit": 2, "insdir": "0 0 1"}
    if explicit_ions:
        molecules = {
            "protein": protein,
            "yttrium": {"atoms": "Yt Cl Cl Cl", "atomic": True, "Ninit": n_yttrium},
            "sodium": {"atoms": "Na", "atomic": True, "Ninit": COUNTERIONS},
        }
    else:
        molecules = {
            "protein": protein,
            "yttrium": {"atoms": "Yt", "atomic": True, "Ninit": n_yttrium},
        }
    translate = {
        name: {"prob": 1.0, "peratom": True} for name in molecules if name != "protein"
    }
    customlj = (
        {
            f"Yt {res}": {"sigma": 2.2, "eps": EPS_Y_CARBOXYLATE}
            for res in ("ASP", "GLU", "CTR")
        }
        if specific
        else {}
    )
    ionic = _ionic_strength(n_yttrium, salt_M, explicit_ions)
    return {
        "customlj": customlj,
        "energy": {"nonbonded": {"coulomb": {"epsr": 78.7, "ionicstrength": ionic}}},
        "atomlist": _atomlist(dp_na, dp_yt),
        "processes": {
            name: {"bound": b, "free": f, "pKd": pka, "pX": PH}
            for name, (b, f, pka) in _TITRATIONS.items()
        },
        "moleculelist": molecules,
        "moves": {
            "atomtranslate": translate,
            "moltransrotcluster": {
                "protein": {
                    "dp": dpp,
                    "dprot": drot,
                    "prob": 1.0,
                    "permol": True,
                    "dir": "0 0 1",
                    "threshold": 6,
                    "clustergroup": "yttrium",
                }
            },
        },
        "system": {
            "temperature": TEMPERATURE,
            "cylinder": {"length": CYLINDER_LENGTH, "radius": CYLINDER_RADIUS},
            "mcloop": {"macro": 10, "micro": micro},
        },
    }
