import numpy as np
import matplotlib.pyplot as plt

# --- configuration ---
M = 201               # grid size in each direction (adjust if needed)
dt = 2.5e-5           # time step
times = [0.0, 0.001, 0.002]

# Transformation for probability plots:
#   None  -> plot p_ij^n
#   "sqrt" -> plot sqrt(p_ij^n)
TRANSFORM = None   # or "sqrt"

def load_field(fname):
    """Load an M×M scalar field saved by Armadillo (raw_ascii)."""
    arr = np.loadtxt(fname)
    # If Armadillo saved as M rows × M cols, this will already be (M,M)
    # but we enforce shape in case of slight differences.
    return arr.reshape((M, M))

def plot_probability(p, t, idx):
    """Plot probability or transformed probability at time t."""
    if TRANSFORM == "sqrt":
        z = np.sqrt(p)
        label = r"$\sqrt{p_{ij}^n}$"
    else:
        z = p
        label = r"$p_{ij}^n$"

    vmax = np.max(z)
    vmin = 0.0

    plt.figure(figsize=(5, 4))
    im = plt.imshow(
        z.T, origin="lower", cmap="viridis",
        extent=[0, 1, 0, 1],
        vmin=vmin, vmax=vmax
    )
    plt.colorbar(im, label=label)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.title(fr"Probability at $t={t:.3f}$ s")
    plt.tight_layout()
    plt.savefig(f"plots/prob_colormap_t{idx}.pdf", dpi=200)
    # plt.show()

def plot_real_imag(re, im, t, idx):
    """Plot Re(u_ij^n) and Im(u_ij^n) at time t."""
    # Use symmetric color scale around zero for both
    max_abs_re = np.max(np.abs(re))
    max_abs_im = np.max(np.abs(im))

    plt.figure(figsize=(5, 4))
    im_re = plt.imshow(
        re.T, origin="lower", cmap="seismic",
        extent=[0, 1, 0, 1],
        vmin=-max_abs_re, vmax= max_abs_re
    )
    plt.colorbar(im_re, label=r"$\Re(u_{ij}^n)$")
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.title(fr"Real part at $t={t:.3f}$ s")
    plt.tight_layout()
    plt.savefig(f"plots/re_colormap_t{idx}.pdf", dpi=200)

    plt.figure(figsize=(5, 4))
    im_im = plt.imshow(
        im.T, origin="lower", cmap="seismic",
        extent=[0, 1, 0, 1],
        vmin=-max_abs_im, vmax= max_abs_im
    )
    plt.colorbar(im_im, label=r"$\Im(u_{ij}^n)$")
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.title(fr"Imaginary part at $t={t:.3f}$ s")
    plt.tight_layout()
    plt.savefig(f"plots/im_colormap_t{idx}.pdf", dpi=200)


def main():
    # t = 0, 0.001, 0.002
    prob_files = [
        "data/prob_t0.dat",
        "data/prob_t1.dat",
        "data/prob_t2.dat",
    ]
    re_files = [
        "data/re_t0.dat",
        "data/re_t1.dat",
        "data/re_t2.dat",
    ]
    im_files = [
        "data/im_t0.dat",
        "data/im_t1.dat",
        "data/im_t2.dat",
    ]

    for idx, (t, f_p, f_r, f_i) in enumerate(zip(times, prob_files, re_files, im_files)):
        p  = load_field(f_p)
        re = load_field(f_r)
        im = load_field(f_i)

        plot_probability(p, t, idx)
        plot_real_imag(re, im, t, idx)

    print("Plots saved in plots/ directory.")

if __name__ == "__main__":
    main()
