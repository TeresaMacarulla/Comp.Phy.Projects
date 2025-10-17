# rz_plot.py
import numpy as np
import matplotlib.pyplot as plt

# Read data from file 
data = np.loadtxt('ry_rx_nointer.txt')
#data = np.loadtxt('ry_rx_yesinter.txt')

# Extract columns
rx1 = data[:, 1]     
ry1 = data[:, 2]
rx2 = data[:, 3]  
ry2 = data[:, 4]   

# Create plot
plt.figure(figsize=(7, 5))
plt.plot(rx1, ry1, color='navy', lw=1.5, label='Particle a')
plt.plot(rx2, ry2, color='darkorange', lw=1.5, label='Particle b')

# Labels and title
plt.xlabel('x-position [μm]', fontsize=12)
plt.ylabel('y-position [μm]', fontsize=12)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

# Save to PDF
plt.tight_layout()
plt.savefig('ryrx_nointer_plot.pdf', format='pdf')

# Optionally display
plt.show()