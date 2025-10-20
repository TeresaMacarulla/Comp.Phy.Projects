import numpy as np
import matplotlib.pyplot as plt

# Initial number of particles
N_initial = 100

# Files and corresponding f values
files = {
    0.1: "trapped_vs_frequency0.100000.txt",
    0.4: "trapped_vs_frequency0.400000.txt",
    0.7: "trapped_vs_frequency0.700000.txt"
}

plt.figure(figsize=(7, 5))

# Loop through each file
for f_val, filename in files.items():
    # Load data (two columns: omegav, n_particles)
    data = np.loadtxt(filename)
    omegav = data[:, 0]
    n_particles = data[:, 1] / N_initial  # normalize

    # Plot normalized number of trapped particles vs frequency
    plt.plot(omegav, n_particles, lw=1.8, label=f"$f = {f_val}$")

# Labels and title
plt.xlabel(r"$\omega_V$ [MHz]", fontsize=12)
plt.ylabel(r"$N_{\mathrm{trapped}} / N_0$", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()

# Save figure as PDF
#plt.savefig("../plots/trapped_vs_frequency.pdf", format='pdf')


# Zoom in frequency------------------------------------------

# Files and corresponding f values
files = {
    0.1: "trapped_vs_freq_zoom0.100000.txt",
    0.4: "trapped_vs_freq_zoom0.400000.txt",
    0.7: "trapped_vs_freq_zoom0.700000.txt"
}

plt.figure(figsize=(7, 5))

# Loop through each file
for f_val, filename in files.items():
    # Load data (two columns: omegav, n_particles)
    data = np.loadtxt(filename)
    omegav = data[:, 0]
    n_particles = data[:, 1] / N_initial  # normalize

    # Plot normalized number of trapped particles vs frequency
    plt.plot(omegav, n_particles, lw=1.8, label=f"$f = {f_val}$")

# Labels and title
plt.xlabel(r"$\omega_V$ [MHz]", fontsize=12)
plt.ylabel(r"$N_{\mathrm{trapped}} / N_0$", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()

# Save figure as PDF
#plt.savefig("../plots/trapped_vs_freq_zoom.pdf", format='pdf')

# Zoom in frequency and interacions------------------------------------------

# Files and corresponding f values
files = {
    0.1: "trap_vs_freq_zoom_inter0.100000.txt",
    0.4: "trap_vs_freq_zoom_inter0.400000.txt",
    0.7: "trap_vs_freq_zoom_inter0.700000.txt"
}

plt.figure(figsize=(7, 5))

# Loop through each file
for f_val, filename in files.items():
    # Load data (two columns: omegav, n_particles)
    data = np.loadtxt(filename)
    omegav = data[:, 0]
    n_particles = data[:, 1] / N_initial  # normalize

    # Plot normalized number of trapped particles vs frequency
    plt.plot(omegav, n_particles, lw=1.8, label=f"$f = {f_val}$")

# Labels and title
plt.xlabel(r"$\omega_V$ [MHz]", fontsize=12)
plt.ylabel(r"$N_{\mathrm{trapped}} / N_0$", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()

# Save figure as PDF
#plt.savefig("../plots/trap_vs_freq_zoom_inter.pdf", format='pdf')

# Optionally show
plt.show()