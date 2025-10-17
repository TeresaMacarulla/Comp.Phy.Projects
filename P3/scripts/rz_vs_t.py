# rz_plot.py
import numpy as np
import matplotlib.pyplot as plt

# Read data from file 
data = np.loadtxt('rz.txt')

# Extract columns
t_us = data[:, 0]     
rz_um = data[:, 1]     

# Create plot
plt.figure(figsize=(7, 5))
plt.plot(t_us, rz_um, color='navy', lw=1.5)

# Labels and title
plt.xlabel('Time [μs]', fontsize=12)
plt.ylabel('z-position [μm]', fontsize=12)
plt.title('Particle motion in Penning Trap', fontsize=13)
plt.grid(True, linestyle='--', alpha=0.6)

# Save to PDF
plt.tight_layout()
# plt.savefig('rz_plot.pdf', format='pdf')

# Optionally display
plt.show()
