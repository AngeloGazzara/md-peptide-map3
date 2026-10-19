"""
Steered MD of the MAP_3 peptide in water: effect of temperature (10 K vs 300 K).

Reads the GROMACS pull outputs (pullf = force, pullx = COM distance) in data/smd/
and plots force and work against the extension.

Run:  python analysis_smd.py
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

DATA = Path(__file__).parent / "data" / "smd"
FIGS = Path(__file__).parent / "figures"
RUNS = {"10 K": "10K", "300 K": "300K"}
COLORS = {"10 K": "#2a78d6", "300 K": "#eb6834"}
KJMOLNM_TO_PN = 1.6605  # 1 kJ mol^-1 nm^-1 = 1.6605 pN


def read_xvg(path):
    """Read a 2-column .xvg file, skipping headers and any incomplete last line."""
    rows = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(("#", "@")):
            continue
        parts = line.split()
        if len(parts) == 2:
            rows.append((float(parts[0]), float(parts[1])))
    return np.array(rows)


def load_run(folder):
    """Return time [ps], pull force [kJ/mol/nm] and COM distance [nm] of one pulling run."""
    time, force = read_xvg(DATA / folder / "pull_pullf.xvg").T
    _, distance = read_xvg(DATA / folder / "pull_pullx.xvg").T
    n = min(len(force), len(distance))
    return time[:n], force[:n], distance[:n]


def main():
    FIGS.mkdir(exist_ok=True)
    fig, (ax_f, ax_w) = plt.subplots(1, 2, figsize=(10, 4))

    for label, folder in RUNS.items():
        time, force, distance = load_run(folder)
        extension = distance - distance[0]                     # ΔL from the initial COM distance
        work = np.concatenate([[0], np.cumsum(0.5 * (force[1:] + force[:-1]) * np.diff(distance))])

        c = COLORS[label]
        ax_f.plot(extension, force, color=c, lw=1, label=label)
        ax_w.plot(extension, work, color=c, lw=2, label=label)
        ax_w.annotate(label, (extension[-1], work[-1]), xytext=(4, 0),
                      textcoords="offset points", va="center", fontsize=9)

        print(f"{label:>6}: {time[-1]:6.0f} ps simulated, ΔL max = {extension.max():.2f} nm, "
              f"F max = {force.max():.0f} kJ/mol/nm ({force.max() * KJMOLNM_TO_PN:.0f} pN), "
              f"total work = {work[-1]:.0f} kJ/mol")

    ax_f.set_xlabel("ΔL [nm]")
    ax_f.set_ylabel("Force [kJ mol$^{-1}$ nm$^{-1}$]")
    ax_f.set_title("Force – extension")
    ax_f.legend(frameon=False)
    ax_w.set_xlabel("ΔL [nm]")
    ax_w.set_ylabel("Work [kJ mol$^{-1}$]")
    ax_w.set_title("Work done by the pulling force")
    for ax in (ax_f, ax_w):
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(alpha=0.25, lw=0.5)
        ax.set_xlim(left=0)

    fig.tight_layout()
    fig.savefig(FIGS / "smd_10K_vs_300K.png", dpi=200)
    print(f"Figure saved to {FIGS / 'smd_10K_vs_300K.png'}")


if __name__ == "__main__":
    main()
