import json
from pathlib import Path

import pytest

from analysis.faunus_input import faunus_input

REFERENCE = Path(__file__).parent / "reference"
LEGACY_STRUCTURE = "../hsa-ph7.aam"  # path the 2017 notebook wrote into hsa.json


def _case(name):
    n, ions, lj, salt = name.removesuffix(".json").split("-")
    return int(n[2:]), float(salt[2:]), ions == "explicit", lj == "specific"


@pytest.mark.parametrize("ref", sorted(REFERENCE.glob("*.json")), ids=lambda p: p.stem)
def test_matches_2017_notebook_input(ref):
    # Reference files were written by the original mkinput() in pmf.ipynb (v1.0.1).
    n, salt, explicit, specific = _case(ref.name)
    got = faunus_input(
        n,
        salt,
        explicit_ions=explicit,
        specific=specific,
        micro=10_000_000,
        structure=LEGACY_STRUCTURE,
    )
    assert got == json.loads(ref.read_text(encoding="utf-8"))


def test_explicit_ions_ignore_added_salt():
    a = faunus_input(16, 0.0, explicit_ions=True, specific=True, micro=1)
    b = faunus_input(16, 0.2, explicit_ions=True, specific=True, micro=1)
    assert a == b


def test_implicit_ionic_strength_grows_with_salt():
    def ionic(salt):
        d = faunus_input(8, salt, explicit_ions=False, specific=True, micro=1)
        return d["energy"]["nonbonded"]["coulomb"]["ionicstrength"]

    assert ionic(0.2) - ionic(0.0) == pytest.approx(0.2)
