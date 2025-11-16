import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.interpolate import make_interp_spline

DATA_DIR = os.path.join("..", "txt")
PLOTS_DIR = os.path.join("..", "plots")
os.makedirs(PLOTS_DIR, exist_ok=True)

sizes = [40, 60, 80, 100]
data = {}

# read files
for L in sizes:
    filename = f"parameters_L_{L}.txt"
    filepath = os.path.join(DATA_DIR, filename)
    arr = np.loadtxt(filepath)

    # Obtain from colummns: T, Cv, chi, eps, mags
    T    = arr[:, 0]
    Cv   = arr[:, 1]
    chi  = arr[:, 2]
    eps  = arr[:, 3]
    mags = arr[:, 4]

    data[L] = {"T": T, "Cv": Cv, "chi": chi, "eps": eps, "mags": mags}

# -----------------------------
# ----- Generate plots ------
# -----------------------------

# 1) Cv vs T
plt.figure()
for L in sizes:
    T  = data[L]["T"]
    Cv = data[L]["Cv"]

    line = plt.plot(T, Cv, linestyle='none', marker='o', markersize=3,
                    label=f"L = {L}")
    color = line[0].get_color()  

    i_max = np.argmax(Cv)
    T_max = T[i_max]
    Cv_max = Cv[i_max]

    plt.plot(T_max, Cv_max, marker='s', markersize=6, color=color)

    plt.annotate(f"({T_max:.2f}, {Cv_max:.3f})",
                 xy=(T_max, Cv_max),
                 xytext=(5, 5), textcoords='offset points',
                 fontsize=8)

plt.xlabel(r"$T$ [J/$k_B$]")
plt.ylabel(r"$C_V/N$")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "Cv_vs_T.pdf"))

# 2) chi vs T
plt.figure()
for L in sizes:
    T   = data[L]["T"]
    chi = data[L]["chi"]

    line = plt.plot(T, chi, linestyle='none', marker='o', markersize=3,
                    label=f"L = {L}")
    color = line[0].get_color()

    i_max = np.argmax(chi)
    T_max = T[i_max]
    chi_max = chi[i_max]

    plt.plot(T_max, chi_max, marker='s', markersize=6, color=color)

    plt.annotate(f"({T_max:.2f}, {chi_max:.3f})",
                 xy=(T_max, chi_max),
                 xytext=(5, 5), textcoords='offset points',
                 fontsize=8)

plt.xlabel(r"$T$ [J/$k_B$]")
plt.ylabel(r"$\chi/N$ [$J^{-1}$]")  
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "chi_vs_T.pdf"))

# 3) eps vs T
plt.figure()
for L in sizes:
    plt.plot(data[L]["T"], data[L]["eps"], linestyle='none', marker='o', markersize=3, label=f"L = {L}")
plt.xlabel(r"$T$ [J/$k_B$]")
plt.ylabel(r"$\langle \epsilon \rangle$ [J]")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "eps_vs_T.pdf"))

# 4) |m| vs T
plt.figure()
for L in sizes:
    plt.plot(data[L]["T"], data[L]["mags"], linestyle='none', marker='o', markersize=3, label=f"L = {L}")
plt.xlabel(r"$T$ [J/$k_B$]")
plt.ylabel(r"$\langle |m| \rangle$")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "mabs_vs_T.pdf"))

