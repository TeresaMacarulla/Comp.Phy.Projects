# rz_plot.py
import numpy as np
import matplotlib.pyplot as plt

# Read data from file 
data4 = np.loadtxt('r_error_4000.txt')
data8 = np.loadtxt('r_error_8000.txt')
data16 = np.loadtxt('r_error_16000.txt')
data32 = np.loadtxt('r_error_32000.txt')

# Extract columns
t4 = data4[:, 0]     
er_RK4_4 = data4[:, 1]    
er_Euler_4 = data4[:, 2] 

t8 = data8[:, 0]     
er_RK4_8 = data8[:, 1]    
er_Euler_8 = data8[:, 2] 

t16 = data16[:, 0]     
er_RK4_16 = data16[:, 1]    
er_Euler_16 = data16[:, 2] 

t32 = data32[:, 0]     
er_RK4_32 = data32[:, 1]    
er_Euler_32 = data32[:, 2] 

# Create plot for RK4
plt.figure(figsize=(7, 5))
plt.plot(t4, er_RK4_4, lw=1.5, label= 'Nsteps=4000')
plt.plot(t8, er_RK4_8, lw=1.5, label= 'Nsteps=8000')
plt.plot(t16, er_RK4_16, lw=1.5, label= 'Nsteps=16000')
plt.plot(t32, er_RK4_32, lw=1.5, label= 'Nsteps=32000')

# Labels and title
plt.xlabel('Time [μs]', fontsize=12)
plt.ylabel('realtive error', fontsize=12)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

# Save to PDF
plt.tight_layout()
plt.savefig('r_error_RK4.pdf', format='pdf')

# Optionally display
plt.show()

# Create plot for Euler
plt.figure(figsize=(7, 5))
plt.plot(t4, er_Euler_4, lw=1.5, label= 'Nsteps=4000')
plt.plot(t8, er_Euler_8, lw=1.5, label= 'Nsteps=8000')
plt.plot(t16, er_Euler_16, lw=1.5, label= 'Nsteps=16000')
plt.plot(t32, er_Euler_32, lw=1.5, label= 'Nsteps=32000')

# Labels and title
plt.xlabel('Time [μs]', fontsize=12)
plt.ylabel('realtive error', fontsize=12)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

# Save to PDF
plt.tight_layout()
plt.savefig('r_error_Euler.pdf', format='pdf')

# Optionally display
plt.show()