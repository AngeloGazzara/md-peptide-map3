"""
Equilibrium MD of MAP_3 in water (2 ns, 300 K): analysis of the GROMACS outputs in data/md/.

Temperature equilibration, backbone RMSD, per-residue RMSF, end-to-end distance,
hydrogen bonds and Ramachandran plot (first vs last frame).

Run:  python analysis_md.py
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from MDAnalysis.analysis.data.filenames import Rama_ref

DATA = Path(__file__).parent / "data" / "md"
FIGS = Path(__file__).parent / "figures"
BLUE, ORANGE, GREY = "#2a78d6", "#eb6834", "#8a8984"


def xvg(name, **kw):
    return np.loadtxt(DATA / name, comments=["#", "@"], unpack=True, **kw)


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(alpha=0.25, lw=0.5)


def main():
    FIGS.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))

    # 1. Temperature: equilibration and distribution over the last ns
    t, temp, _ = xvg("Temp_Press.xvg")
    last_ns = t >= t[-1] - 1000
    ax = axes[0, 0]
    ax.plot(t / 1000, temp, color=GREY, lw=0.6)
    ax.axhline(300, color="k", ls="--", lw=0.8)
    inset = ax.inset_axes([0.55, 0.22, 0.4, 0.36])
    inset.hist(temp[last_ns], bins=25, density=True, color=BLUE, alpha=0.7)
    inset.set_xlabel("T [K] (last ns)", fontsize=8)
    inset.tick_params(labelsize=7)
    inset.set_yticks([])
    style(ax, "Temperature", "Time [ns]", "T [K]")
    print(f"Temperature (last ns): {temp[last_ns].mean():.1f} ± {temp[last_ns].std():.1f} K")

    # 2. Backbone RMSD vs the starting structure and vs the structure at 1.5 ns
    t, rmsd = xvg("rmsd.xvg")
    _, rmsd15 = xvg("rmsd15.xvg")
    ax = axes[0, 1]
    ax.plot(t, rmsd, color=BLUE, lw=0.8, label="ref: 0 ns")
    ax.plot(t, rmsd15, color=ORANGE, lw=0.8, label="ref: 1.5 ns")
    ax.legend(frameon=False)
    ax.set_ylim(bottom=0)
    style(ax, "Backbone RMSD", "Time [ns]", "RMSD [nm]")

    # 3. RMSF per residue (C-alpha)
    res, rmsf = xvg("rmsf.xvg")
    ax = axes[0, 2]
    ax.plot(res, rmsf, color="k", marker="o", ms=4, lw=1)
    ax.set_xticks(res[::2])
    ax.set_ylim(bottom=0)
    style(ax, "C-alpha RMSF", "Residue", "RMSF [nm]")

    # 4. End-to-end distance (COM of residue 1 – COM of residue 18)
    t, d = xvg("distance.xvg")
    ax = axes[1, 0]
    ax.plot(t / 1000, d, color=BLUE, lw=0.6)
    style(ax, "End-to-end distance", "Time [ns]", "$d_{ee}$ [nm]")
    print(f"End-to-end distance: {d.mean():.2f} ± {d.std():.2f} nm")

    # 5. Hydrogen bonds
    t, hb_pw, _ = xvg("hb_water_protein.xvg")
    _, hb_pp, _ = xvg("hb_protein_protein.xvg")
    ax = axes[1, 1]
    ax.plot(t / 1000, hb_pw, color=BLUE, lw=0.6, label="protein–water")
    ax.plot(t / 1000, hb_pp, color=ORANGE, lw=0.6, label="protein–protein")
    ax.legend(frameon=False, loc="center right")
    ax.set_ylim(bottom=0)
    style(ax, "Hydrogen bonds", "Time [ns]", "# H-bonds")
    print(f"H-bonds: protein–water {hb_pw.mean():.1f}, protein–protein {hb_pp.mean():.1f}")

    # 6. Ramachandran plot: first and last frame over the allowed regions
    rows = [line.split() for line in (DATA / "rama.xvg").read_text().splitlines()
            if line.strip() and not line.startswith(("#", "@"))]
    phi = np.array([float(r[0]) for r in rows])
    psi = np.array([float(r[1]) for r in rows])
    residues = np.array([r[2] for r in rows])
    starts = np.flatnonzero(residues == residues[0])      # index where each frame begins
    n = starts[1] - starts[0]                               # residues per frame
    ax = axes[1, 2]
    X, Y = np.meshgrid(np.arange(-180, 180, 4), np.arange(-180, 180, 4))
    ax.contourf(X, Y, np.load(Rama_ref), levels=[1, 17, 15000], colors=["#d6e9fb", "#9cc8f3"])
    for start, color, label in [(starts[0], GREY, "first frame"), (starts[-1], ORANGE, "last frame")]:
        ax.scatter(phi[start:start + n], psi[start:start + n], s=22, color=color,
                   edgecolor="white", lw=0.5, label=label, zorder=3)
    handles = [Line2D([], [], color="#d6e9fb", lw=8), Line2D([], [], color="#9cc8f3", lw=8),
               *ax.get_legend_handles_labels()[0]]
    ax.legend(handles, ["99% allowed", "90% allowed", "first frame", "last frame"],
              frameon=False, fontsize=8, loc="lower right")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(range(-180, 181, 90))
    ax.set_yticks(range(-180, 181, 90))
    style(ax, "Ramachandran plot", "φ [°]", "ψ [°]")

    fig.tight_layout()
    fig.savefig(FIGS / "md_300K.png", dpi=200)
    print(f"Figure saved to {FIGS / 'md_300K.png'}")


if __name__ == "__main__":
    main()
