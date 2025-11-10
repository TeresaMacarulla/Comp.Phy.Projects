
"""
Reads two time-series files with 5 columns each:
1) Monte Carlo cycles
2) epsilon starting from ORDERED initial state
3) running <epsilon> up to that cycle (ORDERED)
4) epsilon starting from UNORDERED (random) initial state
5) running <epsilon> up to that cycle (UNORDERED)

It plots:
- The instantaneous epsilon curves (thin, semi-transparent)
- The running averages (thick), for T=1.0 and T=2.4, ordered and unordered.

Expected files (relative to this script):
  ../txt/epsilon_T_1.000000.txt
  ../txt/epsilon_T_2.400000.txt
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

def load_eps_file(path):
    """
    Returns:
      n_mc, eps_o, mean_eps_o, eps_u, mean_eps_u
    """
    arr = np.loadtxt(path, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 5:
        raise ValueError(f"File {path} must have exactly 5 columns.")
    n_mc = arr[:, 0]
    eps_o = arr[:, 1]
    mean_eps_o = arr[:, 2]
    eps_u = arr[:, 3]
    mean_eps_u = arr[:, 4]
    return n_mc, eps_o, mean_eps_o, eps_u, mean_eps_u

def main():
    script_dir = Path(__file__).resolve().parent
    txt_dir = script_dir.parent / "txt"

    f_T10 = txt_dir / "epsilon_T_1.000000.txt"
    f_T24 = txt_dir / "epsilon_T_2.400000.txt"

    # Load data
    n_MC_cycles_T10, eps_T10_o, mean_eps_T10_o, eps_T10_u, mean_eps_T10_u = load_eps_file(f_T10)
    n_MC_cycles_T24, eps_T24_o, mean_eps_T24_o, eps_T24_u, mean_eps_T24_u = load_eps_file(f_T24)

    # Sanity: reuse a single x-axis if lengths match; otherwise plot with their own x
    if len(n_MC_cycles_T10) != len(n_MC_cycles_T24):
        # Use T=1.0 x-axis for T=1.0 data and T=2.4 x-axis for T=2.4 data
        # For the requested combined plot, we assume same x; if not, we warn but still plot both.
        print("Warning: T=1.0 and T=2.4 have different number of MC cycles; plotting on the same axes anyway.")

    # Build unified x-axis for plotting calls you provided
    # (You asked to use 'n_MC_cycles' directly; we’ll take T=1.0 cycles.)
    n_MC_cycles = n_MC_cycles_T10

    # Create figure
    plt.figure(figsize=(9, 5.5))

    # --- Thin instantaneous epsilon curves (semi-transparent) ---
    plt.plot(n_MC_cycles, eps_T10_u, '-', color='#377eb8', alpha=0.4, linewidth=1.0)
    plt.plot(n_MC_cycles, eps_T10_o, '-', color='#4daf4a', alpha=0.4, linewidth=1.0)
    plt.plot(n_MC_cycles, eps_T24_u, '-', color='#e41a1c', alpha=0.4, linewidth=1.0)
    plt.plot(n_MC_cycles, eps_T24_o, '-', color='#984ea3', alpha=0.4, linewidth=1.0)

    # --- Thick running averages with labels ---
    plt.plot(n_MC_cycles, mean_eps_T10_u, '-', linewidth=2.0, color='#377eb8',
             label=r'<$\epsilon$>, $T=1.0\; J/k_{B}$, unordered')
    plt.plot(n_MC_cycles, mean_eps_T10_o, '-', linewidth=2.0, color='#4daf4a',
             label=r'<$\epsilon$>, $T=1.0\; J/k_{B}$, ordered')
    plt.plot(n_MC_cycles, mean_eps_T24_u, '-', linewidth=2.0, color='#e41a1c',
             label=r'<$\epsilon$>, $T=2.4\; J/k_{B}$, unordered')
    plt.plot(n_MC_cycles, mean_eps_T24_o, '-', linewidth=2.0, color='#984ea3',
             label=r'<$\epsilon$>, $T=2.4\; J/k_{B}$, ordered')

    # Axes, grid, legend
    plt.xlabel("Monte Carlo cycles")
    plt.ylabel(r"Energy per spin $\epsilon$")
    plt.grid(True, alpha=0.25, linestyle='--', linewidth=0.6)
    plt.legend(loc='best', frameon=True)
    plt.xscale('log')
    plt.tight_layout()

    out_pdf = script_dir / "../plots/equilibration_L20.pdf"
    plt.savefig(out_pdf)
    plt.show()

if __name__ == "__main__":
    main()
