import numpy as np
import matplotlib.pyplot as plt
import os

# Directories (relative to this script in 'scripts/')
DATA_DIR = os.path.join("data")
PLOTS_DIR = os.path.join("plots")

# Create plots directory if it does not exist
os.makedirs(PLOTS_DIR, exist_ok=True)

# Full path to the data file
data_file = os.path.join(DATA_DIR, "prob_deviation.txt")

# Load data: first column = time step, second column = probability deviation
data = np.loadtxt(data_file)
time_steps = data[:, 0]
prob_deviation = data[:, 1]

# Create figure
plt.figure()

# Plot deviation vs time step
plt.plot(time_steps, prob_deviation, marker='o', linestyle='none', markersize=3)

plt.xlabel("Time step")
plt.ylabel("Probability deviation from 1")

plt.grid(True, linestyle='--', linewidth=0.5)
plt.tight_layout()

# Save figure as PDF in '../plots'
output_file = os.path.join(PLOTS_DIR, "prob_deviation.pdf")
plt.savefig(output_file)

plt.show()
