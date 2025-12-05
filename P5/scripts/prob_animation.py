import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# ---------------------------------------------------------
# Configuration
# ---------------------------------------------------------
DATA_DIR   = "data/animation"          # folder with prob*.dat
BASE_NAME  = "prob"          # prefix of files
N_FRAMES   = 199             # prob0.dat ... prob198.dat
DT         = 2.5e-5          # time step used in the simulation
OUT_GIF    = "plots/probability_evolution.gif"

# ---------------------------------------------------------
# Collect file names
# ---------------------------------------------------------
file_list = [os.path.join(DATA_DIR, f"{BASE_NAME}{i}.dat") for i in range(N_FRAMES)]

# Load first frame to get grid size
P0 = np.loadtxt(file_list[0])
M_y, M_x = P0.shape   # should both be M

# Compute a global max for a fixed colour scale
global_max = P0.max()
for fname in file_list[1:]:
    P_tmp = np.loadtxt(fname)
    global_max = max(global_max, P_tmp.max())

print(f"Grid size: {M_x} x {M_y}")
print(f"Global max probability: {global_max:.6g}")

# ---------------------------------------------------------
# Set up figure and initial image
# ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(5, 4))

im = ax.imshow(
    P0.T, origin="lower", cmap="viridis",
    extent=[0, 1, 0, 1],
    vmin=0.0, vmax=global_max
)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r"$p_{ij}^n = |u_{ij}^n|^2$")

ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_title(r"$t = 0.0000$")

plt.tight_layout()

# ---------------------------------------------------------
# Animation update function
# ---------------------------------------------------------
def update(frame):
    """Update image for frame index 'frame'."""
    P = np.loadtxt(file_list[frame])
    im.set_data(P.T)

    t = frame * DT
    ax.set_title(rf"$t = {t:.4f}$")

    return [im]

# ---------------------------------------------------------
# Create and save animation as GIF (no ffmpeg needed)
# ---------------------------------------------------------
fps = 30  # frames per second; adjust if needed
anim = FuncAnimation(fig, update, frames=N_FRAMES, blit=True)

os.makedirs("plots", exist_ok=True)

writer = PillowWriter(fps=fps)
anim.save(OUT_GIF, writer=writer)

print(f"Animation saved to {OUT_GIF}")

