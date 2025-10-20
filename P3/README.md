# Penning Trap Simulation Framework

This project implements a numerical simulation of charged particles confined in a Penning Trap.  
It includes simulations of single and multiple particles, with both static and time–dependent potentials, using the Runge–Kutta (RK4) and Forward Euler integration methods.

---

## 🧩 Project Structure

P3/
│
├── include/utils.h/ # Header files (.hpp)
│ ├── Particle.hpp
│ ├── PenningTrap.hpp
│ ├── PenningTrap_Vt.hpp
│ └── Evolution.hpp
│
├── src/ # Source files (.cpp)
│ ├── main.cpp # Main program
│ ├── Particle.cpp # Particle class implementation
│ ├── PenningTrap.cpp # Static trap class (constant V0)
│ ├── PenningTrap_Vt.cpp # Time-dependent trap class
│ └── Evolution.cpp # Auxiliary evolution functions
│
├── scripts/ # Data files and Python plotting scripts
│ ├── *.txt # Simulation outputs
│ ├── *.py # Plotting utilities
│
├── plots/ # Final result figures (.pdf)
│
└── src/main.exe # Compiled executable

---

## ⚙️ Overview

The simulation begins with the construction of a **Penning Trap** through the class `PenningTrap`,  
which defines the physical trap (magnetic field `B0`, electric potential `V0`, and characteristic dimension `d`).  
Each trap contains a vector of `Particle` objects, where each particle is defined by its:
- Charge `q`
- Mass `m`
- Position vector `r`
- Velocity vector `v`

---

## 🧱 Class Descriptions

### **Particle**
Defines a single charged particle.

### **PenningTrap**
Defines the static Penning Trap with constant `B0` and `V0`.  
Contains all the functions required for particle dynamics:

| Function | Description |
|-----------|--------------|
| `add_particle()` | Adds a particle to the trap |
| `external_E_field()` | Returns the electric field vector at position `r` |
| `external_B_field()` | Returns the magnetic field vector at position `r` |
| `force_particle()` | Coulomb force exerted by particle *j* on particle *i* |
| `total_force_external()` | Lorentz force due to external fields |
| `total_force_particles()` | Force on particle *i* from all other particles |
| `total_force()` | Total force on particle *i* (external + interactions) |
| `evolve_RK4()` | RK4 update for all particles over one timestep `h` |
| `evolve_forward_Euler()` | Same as above but using Forward Euler |

---

## 🚀 Workflow — `main.cpp`

The `main.cpp` file orchestrates all simulations in two major parts:

### **Part 1 — Static Potential (constant `V0`)**
1. Runs single–particle and two–particle simulations.
2. Uses methods defined in `PenningTrap.cpp`.
3. Produces data stored as `.txt` files for later plotting in Python.

The following helper functions (from `Evolution.cpp`) are used:

| Function | Purpose |
|-----------|----------|
| `evolution_rz()` | Simulates one particle and outputs `rz.txt` (position vs time) |
| `evolution_two()` | Simulates two particles with/without Coulomb interaction; outputs position (`ry_rx_...`) and velocity (`rx_vx_...`) data |
| `evolution_RK4_Euler()` | Compares RK4 and Euler against the analytical solution; outputs relative errors (`r_error_...txt`) |
| `evolution_nsteps()` | Repeats the simulation for multiple step counts `{4000, 8000, 16000, 32000}` using RK4 and Euler |

All these `.txt` files are later used in Python scripts to generate the plots.

---

### **Part 2 — Time-Dependent Potential \( V(t) = V_0 [1 + f \cos(\omega_v t)] \)**

In this section, the simulation uses the new class **`PenningTrap_Vt`**,  
which extends `PenningTrap` by introducing a time–varying electric potential.

#### Key modifications in `PenningTrap_Vt`:
- The electric field is modulated in time as \( V(t) = V_0 [1 + f \cos(\omega_v t)] \).
- For particles with |r| > d, both **E** and **B** fields are set to zero.
- Adds the function `inside_particles()` that counts how many particles remain inside the trap (|r| < d) at any given time.

---

### **Filling and Evolving the Trap in Part 2**

A helper function `fill_trap()` (defined in `Evolution.cpp`) fills the trap with randomly distributed particles:
- Positions and velocities are drawn from Gaussian distributions.
- The trap typically contains 100 particles.

The trap evolution for each configuration is handled by:
```cpp
evolution_many(trap_Vt, T_us=500.0, Nstep=1000, inter=0);


