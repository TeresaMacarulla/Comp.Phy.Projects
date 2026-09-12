# Computational Physics Projects

This repository contains five numerical physics projects covering a broad range of computational methods, from finite-difference solvers and eigenvalue algorithms to particle simulations, Monte Carlo methods, and time-dependent quantum mechanics.

The projects are organized independently in the `P1/`–`P5/` folders. Each folder contains its own README with more detailed information about the numerical method, implementation, workflow, data files, and plotting scripts.

## Repository structure

```text
Comp.Phy.Projects/
├── P1/
│   └── README.md
├── P2/
│   └── README.md
├── P3/
│   └── README.md
├── P4/
│   └── README.md
├── P5/
│   └── README.md
└── README.md
```

## Project 1 — Finite-difference solution of a boundary-value problem

`P1/` introduces numerical solution methods for a one-dimensional boundary-value problem.

The project includes:

- generation of the exact solution for comparison,
- implementation of a general finite-difference algorithm,
- implementation of a specialized algorithm for the same tridiagonal system,
- absolute and relative error analysis,
- comparison of execution times for the general and specialized methods,
- C++ programs for the numerical calculations,
- Python scripts for visualization and post-processing.

The project studies how the numerical approximation changes with increasing discretization size and compares accuracy and computational efficiency between the two algorithms.

For detailed file descriptions and workflow, see:

**[`P1/README.md`](P1/README.md)**

## Project 2 — Eigenvalue problems and the Jacobi rotation method

`P2/` focuses on numerical methods for symmetric eigenvalue problems.

The project includes:

- construction of a tridiagonal matrix,
- comparison between analytical and numerical eigenvalues and eigenvectors,
- implementation of a routine for locating the largest off-diagonal matrix element,
- implementation of the **Jacobi rotation algorithm**,
- investigation of the number of transformations required as matrix size increases,
- calculation of the three lowest eigenstates,
- comparison between numerical and analytical eigenvectors,
- Python scripts for plotting and analysis.

The implementation uses **C++** for the numerical algorithms and **Python** for visualization.

For the complete workflow and file descriptions, see:

**[`P2/README.md`](P2/README.md)**

## Project 3 — Charged particles in a Penning trap

`P3/` implements a numerical framework for simulating charged particles confined in a **Penning trap**.

The project models both single- and multi-particle systems and includes:

- electric and magnetic fields,
- Coulomb interactions between particles,
- **Runge–Kutta 4 (RK4)** and **Forward Euler** time integration,
- comparison with analytical solutions,
- convergence and relative-error studies,
- static and time-dependent trapping potentials,
- simulations of many particles with randomized initial conditions,
- analysis of how many particles remain confined in the trap.

The code is organized around `Particle`, `PenningTrap`, and time-dependent `PenningTrap_Vt` classes, with additional evolution functions for the different numerical experiments.

For detailed information about the classes, simulation workflow, and generated data, see:

**[`P3/README.md`](P3/README.md)**

## Project 4 — 2D Ising model with Monte Carlo methods

`P4/` studies the two-dimensional **Ising model** using the Metropolis Markov Chain Monte Carlo algorithm.

The project includes:

- initialization of ordered and disordered spin lattices,
- Metropolis sampling with periodic boundary conditions,
- calculation of energy and magnetization,
- heat capacity and magnetic susceptibility,
- equilibration and burn-in studies,
- energy-distribution histograms,
- temperature scans for different lattice sizes,
- **OpenMP parallelization** across independent temperatures,
- finite-size scaling and extrapolation of the critical temperature.

The numerical simulations are implemented in **C++**, while **Python**, NumPy, and Matplotlib are used for statistical analysis and plotting.

For build instructions, simulation options, data formats, and analysis scripts, see:

**[`P4/README.md`](P4/README.md)**

## Project 5 — 2D time-dependent Schrödinger equation

`P5/` implements a numerical solver for the two-dimensional, time-dependent **Schrödinger equation** using the **Crank–Nicolson method**.

The solver is used to study diffraction and interference of a Gaussian wave packet interacting with one-, two-, and three-slit potentials.

The project includes:

- sparse-matrix construction for the Crank–Nicolson scheme,
- initialization of slit potentials,
- normalized Gaussian wave packets,
- time evolution of the quantum state,
- probability-conservation analysis,
- visualization of probability density and the real and imaginary parts of the wave function,
- screen-detection probability and interference patterns,
- animated probability-density evolution.

The numerical solver is written in **C++** using Armadillo, while **Python** is used for post-processing, figures, and animations.

For the full interactive workflow, numerical functions, and visualization tools, see:

**[`P5/README.md`](P5/README.md)**

## Computational methods covered

Across the five projects, the repository explores several core techniques in computational physics:

- Finite-difference discretization
- Tridiagonal linear systems
- Error and convergence analysis
- Jacobi rotations for eigenvalue problems
- Numerical eigenvectors and eigenstates
- Forward Euler integration
- Fourth-order Runge–Kutta integration
- Classical charged-particle dynamics
- Metropolis Markov Chain Monte Carlo
- Statistical-mechanics observables
- Parallel computing with OpenMP
- Finite-size scaling
- Sparse linear algebra
- Crank–Nicolson time evolution
- Numerical solutions of the Schrödinger equation

## Main technologies

- **C++**
- **Python**
- **Armadillo**
- **OpenMP**
- **NumPy**
- **Matplotlib**

## Documentation

This root README provides a high-level overview of the repository. More detailed documentation is available inside each project folder:

- **[P1 README](P1/README.md)** — finite-difference algorithms, error analysis, and performance comparison.
- **[P2 README](P2/README.md)** — eigenvalue problems and the Jacobi rotation method.
- **[P3 README](P3/README.md)** — Penning-trap particle simulations and numerical time integration.
- **[P4 README](P4/README.md)** — 2D Ising model, Metropolis MCMC, OpenMP, and finite-size scaling.
- **[P5 README](P5/README.md)** — Crank–Nicolson solution of the 2D time-dependent Schrödinger equation and slit-interference simulations.
