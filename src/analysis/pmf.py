"""Protein-protein potential of mean force from the two-body Monte Carlo runs."""

import csv
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import UnivariateSpline

AVOGADRO = 6.02214076e23  # 1/mol
CYLINDER_LENGTH = 563.0  # Å
CYLINDER_RADIUS = 50.0  # Å
CYLINDER_VOLUME = math.pi * CYLINDER_RADIUS**2 * CYLINDER_LENGTH  # Å³

# Spline smoothing per Y³⁺ count, hand-tuned in 2017; kept so figures match the paper.
SMOOTHING = {1: 40, 8: 40, 16: 60, 160: 40}

Curves = dict[tuple[float, int], tuple[np.ndarray, np.ndarray]]


def pmf_from_rdf(r, g, r_tail=250.0):
    """Potential of mean force, w(r)/kT = -ln g(r).

    Parameters
    ----------
    r : array_like
        Mass-centre separation (Å).
    g : array_like
        Unnormalized distance histogram along the cylinder axis.
    r_tail : float
        Beyond this separation (Å) the proteins no longer interact, so the
        histogram is scaled to average one there, which sets w(∞) = 0.

    Returns
    -------
    numpy.ndarray
        w(r) in units of kT.
    """
    r, g = np.asarray(r), np.asarray(g)
    return -np.log(g / g[r > r_tail].mean())


def yttrium_mM(n_yttrium):
    """Y³⁺ concentration (mM) for `n_yttrium` ions in the simulation cylinder."""
    litres = CYLINDER_VOLUME * 1e-27
    return n_yttrium / (AVOGADRO * litres) * 1e3


def load_pmf(path):
    """Read PMFs from the long-format CSV in ``data/processed``.

    Parameters
    ----------
    path : str or pathlib.Path
        CSV with columns ``salt_M, n_yttrium, r_angstrom, pmf_kT``.

    Returns
    -------
    dict
        ``(salt_M, n_yttrium) -> (r, w)`` with r in Å and w in kT, sorted by r.
    """
    data = np.genfromtxt(Path(path), delimiter=",", names=True, encoding="utf-8")
    curves: Curves = {}
    for salt, n in sorted(
        {(float(s), int(k)) for s, k in data[["salt_M", "n_yttrium"]]}
    ):
        rows = data[(data["salt_M"] == salt) & (data["n_yttrium"] == n)]
        rows = rows[np.argsort(rows["r_angstrom"])]
        curves[(salt, n)] = (rows["r_angstrom"], rows["pmf_kT"])
    return curves


def write_pmf(path, curves):
    """Write PMFs in the long CSV format read by `load_pmf`.

    Parameters
    ----------
    path : str or pathlib.Path
        Output file.
    curves : dict
        ``(salt_M, n_yttrium) -> (r, w)``, as returned by `load_pmf`.
    """
    with Path(path).open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["salt_M", "n_yttrium", "r_angstrom", "pmf_kT"])
        for (salt, n), (r, w) in sorted(curves.items()):
            # repr keeps full float precision, so the file round-trips exactly
            writer.writerows(
                [float(salt), int(n), repr(float(a)), repr(float(b))]
                for a, b in zip(r, w, strict=True)
            )


def smooth_pmf(r, w, n_yttrium):
    """Cubic smoothing spline of a PMF, using the published smoothing factor.

    Returns
    -------
    scipy.interpolate.UnivariateSpline
        Callable w(r).
    """
    spline = UnivariateSpline(r, w, k=3)
    spline.set_smoothing_factor(SMOOTHING[n_yttrium])
    return spline


def plot_pmf(ax, curves, salt_M):
    """Plot raw (faint) and smoothed PMFs for one added-salt concentration.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    curves : dict
        Output of `load_pmf`.
    salt_M : float
        Added 1:1 salt (mol/l) to select.
    """
    colors = ["black", "red", "green", "blue"]
    r_smooth = np.linspace(45, 160, 1000)
    selected = sorted((n, rw) for (s, n), rw in curves.items() if s == salt_M)
    for color, (n, (r, w)) in zip(colors, selected, strict=True):
        label = f"$c_s$={yttrium_mM(n):.0f} mM"
        ax.plot(r, w, "-", color=color, lw=2.5, zorder=1, alpha=0.4)
        ax.plot(
            r_smooth,
            smooth_pmf(r, w, n)(r_smooth),
            color=color,
            lw=2,
            zorder=2,
            label=label,
        )
    ax.set_xlabel(r"$r$/Å")
    ax.set_ylabel(r"$w(r)/k_BT$")
    ax.set_xlim(35, 160)
    ax.set_ylim(-3, 4)
    ax.legend(loc=0, frameon=False)
