from pathlib import Path

import numpy as np
import pytest

from analysis import pmf

ROOT = Path(__file__).parents[1]
CSV = ROOT / "data" / "processed" / "pmf.csv"


def test_pmf_from_rdf_recovers_potential_and_ignores_tail_scale():
    r = np.linspace(40, 300, 1301)
    w_true = -2.0 * np.exp(-(((r - 60) / 10) ** 2))
    g = 3.7 * np.exp(-w_true)  # arbitrary histogram scale
    np.testing.assert_allclose(pmf.pmf_from_rdf(r, g), w_true, atol=1e-12)


def test_yttrium_concentration_matches_published_labels():
    labels = [f"{pmf.yttrium_mM(n):.0f}" for n in (1, 8, 16, 160)]
    assert labels == ["0", "3", "6", "60"]
    assert pmf.yttrium_mM(160) == pytest.approx(60.08, abs=0.01)


# Characterization: pinned from the 2017 pickle (monte-carlo-data.p, Zenodo v1.0.1).
LEGACY = {
    (0.0, 1): (1207, -1.3100915695761544, [-0.737759, 0.390986, 0.318776]),
    (0.0, 8): (1204, -2.6447807958181553, [-1.510888, -0.501525, 0.019439]),
    (0.0, 16): (1203, -2.186683726947077, [-0.422301, 0.001329, 0.077927]),
    (0.0, 160): (1213, -0.7030050558059409, [1.136393, -0.008835, -0.015475]),
    (0.1, 8): (1207, -1.8938481820194273, [-1.253503, -0.34651, 0.028215]),
    (0.2, 1): (1207, -1.2839578494678456, [-0.698152, -0.25085, 0.008368]),
    (0.2, 16): (1209, -1.8066003209157566, [-1.219555, -0.368777, 0.021933]),
    (0.2, 160): (1209, -0.9449685071255506, [-0.332771, -0.273441, 0.024512]),
}


def test_load_pmf_matches_legacy_data():
    curves = pmf.load_pmf(CSV)
    assert len(curves) == 12
    for key, (npoints, wmin, _) in LEGACY.items():
        r, w = curves[key]
        assert len(r) == npoints
        assert w.min() == wmin
        assert np.all(np.diff(r) > 0)


def test_smoothed_pmf_matches_legacy_spline():
    curves = pmf.load_pmf(CSV)
    for (salt, n), (_, _, expected) in LEGACY.items():
        r, w = curves[(salt, n)]
        spline = pmf.smooth_pmf(r, w, n)
        np.testing.assert_allclose(spline([50.0, 80.0, 120.0]), expected, atol=1e-6)


def test_plot_pmf_labels_by_yttrium_concentration():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    pmf.plot_pmf(ax, pmf.load_pmf(CSV), salt_M=0.2)
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels == [f"$c_s$={c} mM" for c in (0, 3, 6, 60)]
    plt.close(fig)


def test_write_pmf_round_trips(tmp_path):
    curves = pmf.load_pmf(CSV)
    out = tmp_path / "pmf.csv"
    pmf.write_pmf(out, curves)
    assert out.read_text(encoding="utf-8") == CSV.read_text(encoding="utf-8")
