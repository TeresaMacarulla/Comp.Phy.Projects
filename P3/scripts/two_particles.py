# rz_plot.py
import numpy as np
import matplotlib.pyplot as plt

# Read data from files 
data_no = np.loadtxt('ry_rx_nointer.txt')
data_yes = np.loadtxt('ry_rx_yesinter.txt')

data_no_vx = np.loadtxt('rx_vx_nointer.txt')
data_yes_vx = np.loadtxt('rx_vx_yesinter.txt')

data_no_vz = np.loadtxt('rz_vz_nointer.txt')
data_yes_vz = np.loadtxt('rz_vz_yesinter.txt')

# Extract columns
rx1_no = data_no[:, 1]     
ry1_no = data_no[:, 2]
rx2_no = data_no[:, 3]  
ry2_no = data_no[:, 4]  

rx1_yes = data_yes[:, 1]     
ry1_yes = data_yes[:, 2]
rx2_yes = data_yes[:, 3]  
ry2_yes = data_yes[:, 4]  
   
vx1_no = data_no_vx[:, 2]
vx2_no = data_no_vx[:, 4]  

rz1_no = data_no_vz[:, 1]     
vz1_no = data_no_vz[:, 2]
rz2_no = data_no_vz[:, 3]  
vz2_no = data_no_vz[:, 4]  

vx1_yes = data_yes_vx[:, 2]
vx2_yes = data_yes_vx[:, 4]  

rz1_yes = data_yes_vz[:, 1]     
vz1_yes = data_yes_vz[:, 2]
rz2_yes = data_yes_vz[:, 3]  
vz2_yes = data_yes_vz[:, 4]  

# Create plot ry vs rx for no interaction -----------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(rx1_no, ry1_no, color='navy', lw=1.5, label='Particle a')
plt.plot(rx2_no, ry2_no, color='darkorange', lw=1.5, label='Particle b')

# Labels and title
plt.xlabel('x-position [μm]', fontsize=12)
plt.ylabel('y-position [μm]', fontsize=12)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

# Save to PDF
plt.tight_layout()
#plt.savefig('../plots/ryrx_nointer_plot.pdf', format='pdf')


# Create plot ry vs rx for interaction -----------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(rx1_yes, ry1_yes, color='navy', lw=1.5, label='Particle a')
plt.plot(rx2_yes, ry2_yes, color='darkorange', lw=1.5, label='Particle b')

# Labels and title
plt.xlabel('x-position [μm]', fontsize=12)
plt.ylabel('y-position [μm]', fontsize=12)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

# Save to PDF
plt.tight_layout()
#plt.savefig('../plots/ryrx_yesinter_plot.pdf', format='pdf')


# Create plot vx vs rx for no interaction -----------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(rx1_no, vx1_no, color='navy', lw=1.5, label='Particle a')
plt.plot(rx2_no, vx2_no, color='darkorange', lw=1.5, label='Particle b')

# Labels and title
plt.xlabel('x-position [μm]', fontsize=12)
plt.ylabel('vx-velocity [μm/μs]', fontsize=12)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

# Save to PDF
plt.tight_layout()
#plt.savefig('../plots/vxrx_nointer_plot.pdf', format='pdf')


# Create plot vz vs rz for no interaction -----------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(rz1_no, vz1_no, color='navy', lw=1.5, label='Particle a')
plt.plot(rz2_no, vz2_no, color='darkorange', lw=1.5, label='Particle b')

# Labels and title
plt.xlabel('z-position [μm]', fontsize=12)
plt.ylabel('vz-velocity [μm/μs]', fontsize=12)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

# Save to PDF
plt.tight_layout()
#plt.savefig('../plots/vzrz_nointer_plot.pdf', format='pdf')


# Create plot vx vs rx for interaction -----------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(rx1_yes, vx1_yes, color='navy', lw=1.5, label='Particle a')
plt.plot(rx2_yes, vx2_yes, color='darkorange', lw=1.5, label='Particle b')

# Labels and title
plt.xlabel('x-position [μm]', fontsize=12)
plt.ylabel('vx-velocity [μm/μs]', fontsize=12)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

# Save to PDF
plt.tight_layout()
#plt.savefig('../plots/vxrx_yesinter_plot.pdf', format='pdf')


# Create plot vz vs rz for interaction -----------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(rz1_yes, vz1_yes, color='navy', lw=1.5, label='Particle a')
plt.plot(rz2_yes, vz2_yes, color='darkorange', lw=1.5, label='Particle b')

# Labels and title
plt.xlabel('z-position [μm]', fontsize=12)
plt.ylabel('vz-velocity [μm/μs]', fontsize=12)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.6)

# Save to PDF
plt.tight_layout()
#plt.savefig('../plots/vzrz_yesinter_plot.pdf', format='pdf')

# Optionally display
plt.show()