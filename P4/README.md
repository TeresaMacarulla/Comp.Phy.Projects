# 2D Ising Model – Monte Carlo & Finite-Size Scaling

C++ + OpenMP implementation of the 2D Ising model with Metropolis MCMC, plus Python
scripts for analysis and plotting.  
The code reproduces all simulations and figures used in the report.

---

## 1. Project Structure

```text
P4/
├── include/
│   └── functions.hpp          # Declarations for the core C++ functions
├── src/
│   ├── main.cpp               # Main program + menu of simulations
│   └── utils/
│       └── functions.cpp      # Implementation of Ising model routines
├── scripts/                   # Python analysis & plotting tools
│   ├── equilibration_results.py
│   ├── fit_Tc_vs_L.py
│   ├── histogram.py
│   └── plot_parameters.py
├── txt/                       # Output data from C++ runs (.txt)
│   ├── eps_T_*.txt
│   ├── eps_T_*_burn-in_*.txt
│   └── parameters_L_*.txt
└── plots/                     # Figures produced by Python scripts (.pdf)
    ├── equilibration_L20.pdf
    ├── hist_eps_T_*.pdf
    ├── eps_vs_T.pdf
    ├── mabs_vs_T.pdf
    ├── Cv_vs_T.pdf
    ├── chi_vs_T.pdf
    └── Tc_vs_invL_fit.pdf
````

---

## 2. Core C++ Functions (`src/utils/functions.cpp`)

These routines implement the Ising model physics and the Monte Carlo algorithm:

* **`initialize_lattice(L, ordered=false)`**
  Creates an `L × L` spin lattice:

  * `ordered = true`: all spins aligned (`+1`).
  * `ordered = false`: random spins (`±1` with equal probability).

* **`compute_energy_and_magnetization(lattice, L, J)`**
  Computes **E** (total energy) and **M** (total magnetisation) for a given microstate,
  using periodic boundary conditions.

* **`run_mcmc_metropolis(lattice, L, T, J, MC_sweeps, burn_in, ...)`**
  Runs the **Metropolis Markov Chain Monte Carlo**:

  * Single–spin flips with precomputed Boltzmann factors.
  * Periodic boundary conditions.
  * Discards the first `burn_in` sweeps.
  * Returns the final lattice and time series of **E** and **M**.

* **`compute_parameters(results, L, T)`**
  From the MCMC samples (typically the output of `run_mcmc_metropolis`), computes:

  * Mean energy per spin ⟨ε⟩,
  * Mean absolute magnetisation per spin ⟨|m|⟩,
  * Heat capacity per spin (C_V/N),
  * Susceptibility per spin (\chi/N).

These functions are declared in `include/functions.hpp` and used by `src/main.cpp`.

---

## 3. Main Program and Simulation Menu (`src/main.cpp`)

Compiling and running `main.cpp` launches an interactive menu in the terminal:

```text
=== ISING MODEL SIMULATION MENU ===
 Choose an option
  1. Apply the Markov Chain Monte Carlo approach to sample spin
     configurations and compute mean energy and magnetisation,
     heat capacity and susceptibility per spin for a fixed temperature.
  2. Study the burn-in time (or equilibration time).
  3. Run a simulation choosing the adequate burn-in Monte Carlo cycles.
  4. Choose some temperature values and run simulations with and without
     parallelization (OpenMP). Compare the timing.
  5. Run parallelization from T = 2.1 to 2.4 J/kB and from L = 40 to 100.
```

Roughly, the options correspond to the report’s Problems:

1. **Problem 4 / 2×2 & fixed-T observables**

   * Uses `run_mcmc_metropolis` + `compute_parameters` for a single temperature.
   * Can be used to validate against analytical results (e.g. (L=2)).

2. **Problem 5 / Equilibration (burn-in) for L=20**

   * Records energy per spin vs. Monte Carlo cycle.
   * Produces `.txt` files later used by `equilibration_results.py`.

3. **Problem 6 / Energy distribution**

   * Runs at a fixed temperature with a chosen burn-in.
   * Saves ε samples to `.txt` for histogram analysis.

4. **Parallel vs Serial timing (OpenMP)**

   * Runs a set of user-chosen temperatures both:

     * in **parallel** (OpenMP `#pragma omp parallel for` over temperatures), and
     * in **serial** (plain `for` loop).
   * Prints timing and speed-up factor.

5. **Full parallel scan in T and L (finite-size scaling)**

   * Uses OpenMP to sweep temperatures from **2.1 to 2.4 J/kB**.
   * For each (T, L): 100 000 Monte Carlo sweeps, 1 000 burn-in sweeps.
   * Stores per-temperature observables to `txt/parameters_L_*.txt`.

Running the options **in order** with appropriate choices of parameters reproduces all
numerical results used in the report.

---

## 4. Data Files (`txt/`)

C++ simulations write their output as plain text files:

* **Equilibration runs**

  * `epsilon_T_1.000000.txt`
  * `epssion_T_2.400000.txt`
    Each file:
    `MC_cycle   epsilon_ordered   <epsilon>_ordered   epsilon_unordered   <epsilon>_unordered`

* **Energy histograms**

  * `eps_T_1.000000_burn-in_100.txt`, `epsilon_T_2.400000_burn-in_100000.txt`, …
    Two-column files (MC cycle, ε) used to build (p_\epsilon(\epsilon; T)).

* **Thermodynamic observables vs T**

  * `parameters_L_40.txt`, `parameters_L_60.txt`, `parameters_L_80.txt`, `parameters_L_100.txt`
    Columns typically contain:
    `T   Cv/N   chi/N   <epsilon>   <|m|>  ` (per lattice size).

These files are the input for the Python plotting scripts.

---

## 5. Python Analysis & Plotting (`scripts/` → `plots/`)

The Python scripts use `NumPy` and `Matplotlib` to turn `.txt` data into publication-ready figures.

* **`equilibration_results.py`**

  * Reads `epsilon_T_*.txt`.
  * Plots instantaneous ε and running ⟨ε⟩ vs Monte Carlo cycles for ordered/disordered starts.
  * Output: `plots/equilibration_L20.pdf`.

* **`histogram.py`**

  * Builds **normalised histograms** of ε (probability per discrete ε level).
  * Output: `plots/hist_eps_T_*.pdf`.

* **`plot_parameters.py`**

  * Reads `parameters_L_*.txt`.
  * Produces:

    * `eps_vs_T.pdf` (energy per spin),
    * `mabs_vs_T.pdf` (absolute magnetisation per spin),
    * `Cv_vs_T.pdf` (heat capacity per spin),
    * `chi_vs_T.pdf` (susceptibility per spin),
      for all lattice sizes.

* **`fit_Tc_vs_L.py`**

  * Uses the pseudo-critical temperatures (T_c(L)) from Table III in the report.
  * Performs a weighted linear fit of (T_c(L)) vs (1/L):
    (T_c(L) = a/L + T_c(\infty)).
  * Output: `plots/Tc_vs_invL_fit.pdf` with fit line and fitted parameters.

All figures included in the report (equilibration, histograms, observables vs T, and
finite-size extrapolation of (T_c)) are generated by these scripts from the data
produced by the C++ simulations.

---

## 6. Typical Workflow

1. **Compile the C++ code** (with OpenMP enabled), e.g.

   ```bash
   g++ -std=c++17 -O2 -I include src/main.cpp src/utils/*.cpp -larmadillo -o src/main.exe -fopenmp
   ```

2. **Run `main.exe`** and follow the menu:

   1. Validate MCMC and observables at fixed T (e.g. (L=2)).
   2. Study burn-in for (L=20) and T = 1.0, 2.4 .
   3. Generate long ε time series for histograms.
   4. Test OpenMP speed-up for selected temperatures at fixed L.
   5. Run full parallel scans in T and L = 40–100 → generates `parameters_L_*.txt`.

3. **Post-process with Python**:

   ```bash
   cd P4/scripts
   python equilibration_results.py
   python histogram.py
   python plot_parameters.py
   python fit_Tc_vs_L.py
   ```

   (Figures will be saved to `P4/plots`.)

4. **Use the PDFs in `plots/`** directly in the LaTeX report.

---

## 7. Notes

* All simulations use units with (k_B = 1) and coupling (J = 1), so temperature is reported in (J/k_B).
* OpenMP parallelisation is applied only over **independent temperatures**.

```
