

## P5 – 2D Time-Dependent Schrödinger Equation (Crank–Nicolson Double Slit)

This project implements a numerical solver for the two-dimensional, time-dependent Schrödinger equation for a single particle in a box with Dirichlet boundary conditions. The equation is discretised with a Crank–Nicolson scheme and used to study diffraction and interference of a Gaussian wave packet on one-, two-, and three-slit potentials.

The code is written in C++ (time integration and linear algebra) and Python (post-processing and plotting).

---

## 1. Folder structure

```text
P5/
  data/                   # numerical output (.dat, .txt files)
    animation/            # (optional) per-frame probability snapshots for video
    im_t0.dat, im_t1.dat, im_t2.dat
    prob_t0.dat, prob_t1.dat, prob_t2.dat
    prob_deviation.txt, prob_deviation_7*.txt
    re_t0.dat, re_t1.dat, re_t2.dat
    screen_prob_t2.dat
    ...
  include/
    functions.hpp         # C++ function declarations
  plots/                  # figures, PDFs, GIFs, animations (created by Python)
  scripts/
    plot_wavefunction_fields.py  # colormaps for probability, Re(u), Im(u)
    prob_animation.py            # animation of probability evolution
    prob_deviation.py            # probability deviation vs time
    screen_probability.py        # detection probability on a screen
  src/
    main.cpp              # main program (interactive, menu-driven)
    utils/
      functions.cpp       # numerical kernels, matrix assembly, helpers
  README.md               # this file
````

---

## 2. Building the C++ code

From inside `P5/`:

```bash
g++ -std=c++17 -O2 -I include src/main.cpp src/utils/*.cpp -larmadillo -o src/main.exe -fopenmp
```

---

## 3. main.cpp – interactive, menu-driven workflow

Run the program from `P5/`:

```bash
src/main.exe
```

The **interactive nature of `main.cpp`** means you do not have to remember function calls or file names: the terminal guides you through the simulation setup and tells you which Python script to run to visualise the results.

### 3.1. Top-level menu

After starting, you see the first menu:

```text
=== Double-slit setup ===
 Choose an option (write a number)
  1. Create A and B matrices and check.
  2. Initialise the potential.
  3. Run a simulation.
```

* **Option 1 – Create A and B matrices and check**
  For debugging and verification. This constructs the Crank–Nicolson matrices (A) and (B) for a given grid size and prints their sparsity patterns using `print_sp_matrix_structure`. Helpful to confirm the stencil structure.

* **Option 2 – Initialise the potential**
  Lets you test and visualise the potential (V(x,y)) (e.g. single, double, triple slit) without running a full time evolution. You provide the potential parameters and can:

  * Print a coarse ASCII representation to the terminal, and
  * Save a PGM image of the potential.

* **Option 3 – Run a simulation**  🚀
  This is the main, most interesting option. It:

  1. Asks for simulation parameters (grid, time step, total time),
  2. Asks for potential characteristics (wall thickness, position, number of slits, etc.),
  3. Asks for initial Gaussian wave packet parameters ((x_c, y_c, \sigma_x, \sigma_y, p_x, p_y)),
  4. Builds the matrices and initial state,
  5. Enters a time-stepping loop.

  After these inputs, you see a **second menu** that selects what kind of analysis/output you want from the simulation.

### 3.2. Simulation submenu (where the “magic” happens)

Once you choose *“Run a simulation”* and the setup is complete, you get:

```text
Now what? Choose an option
 1. Study the deviation of the total probability from 1.0.
 2. Study time evolution of the 2D probability function.
 3. Analyze the detections on a screen.
 4. Create an animation of your simulation.
```

These four options control how data is collected and which `.dat` or `.txt` files are written. The terminal will also remind you which Python script to run afterwards.

#### Option 1 – Probability conservation

* Computes the total probability at each time step, and stores the deviation (|P^n - 1|) in `data/prob_deviation.txt`.
* After the run the program prints a hint:

  ```text
  Run the python script from scripts/prob_deviation.py to see your probability deviations results
  ```
* `scripts/prob_deviation.py` then plots the deviations versus time and saves the figure under `plots/`.

#### Option 2 – 2D probability and wavefunction snapshots

* Stores the full 2D fields at selected times (typically (t=0), (t=0.001), (t=0.002)):

  * Probability: `prob_t0.dat`, `prob_t1.dat`, `prob_t2.dat`
  * Real part: `re_t0.dat`, `re_t1.dat`, `re_t2.dat`
  * Imaginary part: `im_t0.dat`, `im_t1.dat`, `im_t2.dat`
* At the end, the program prints a reminder:

  ```text
  You can find your .dat files in 'data'. Run scripts/plot_wavefunction_fields.py
  to generate the colourmaps of the probability, Re(u_ij) and Im(u_ij).
  ```
* `scripts/plot_wavefunction_fields.py` reads these files and produces colormaps saved in `plots/` (PDF), showing the evolution of the wave packet.

#### Option 3 – Screen detection probability (Problem 9)

* Mimics a detector screen at (x = 0.8) at time (t = 0.002).
* Extracts the vertical slice (u(x=0.8,y,t=0.002)), computes the probability and normalises it so that the 1D distribution sums to 1.
* Stores the result in `data/screen_probability_t2.dat` with two columns:

  ```text
  y   p(y | x=0.8, t=0.002)
  ```
* The terminal then instructs you to run:

  ```text
  Run scripts/screen_probability.py to visualise the detection probability.
  ```
* `scripts/screen_probability.py` produces a 1D plot of the interference fringes on the detection screen.

#### Option 4 – Animation of the probability evolution

* At each time step, saves the 2D probability field into a series of files (e.g. `prob0.dat`, `prob1.dat`, …) in `data/animation` folder.
* After the run, the program points you to:

  ```text
  Run scripts/prob_animation.py to create a GIF of the simulation.
  ```
* `scripts/prob_animation.py` builds an animated GIF (`plots/probability_evolution.gif`) that shows the wave packet moving, diffracting, and interfering over time.

---

## 4. Typical workflow

A common sequence to reproduce the double-slit experiments from the report:

1. **Build the code**

2. **Start the program**

3. **At the first menu**, choose `3` (Run a simulation).

4. **Provide numerical and physical parameters** when prompted:

   * Grid dimension `M` (e.g. 201),
   * Time step `dt` (e.g. `2.5e-5`),
   * Total simulation time `T` (e.g. `0.002`),
   * Potential parameters (`v0`, wall thickness, position, slit aperture, separation, number of slits),
   * Initial Gaussian packet parameters (`xc`, `yc`, `sigma_x`, `sigma_y`, `p_x`, `p_y`).

5. **At the simulation submenu**, select the analysis you want:

   * `1` → generate `prob_deviation.txt` and then run `python scripts/prob_deviation.py`;
   * `2` → generate `prob_t*.dat`, `re_t*.dat`, `im_t*.dat` and run `python scripts/plot_wavefunction_fields.py`;
   * `3` → generate `screen_prob_t2.dat` and run `python scripts/screen_probability.py`;
   * `4` → generate a sequence of `prob*.dat` files and run `python scripts/prob_animation.py`.

The interactive prompts in the terminal give you all the information you need at each step and explicitly tell you which Python script to run to visualise the results.

---

## 5. Summary of key C++ functions

A deeper explanation of each function is given in `src/utils/functions.cpp`.
Here is a brief overview, organised roughly in the order they enter the workflow of `main.cpp`.

* **Map internal grid indices ((i,j)) to a single vector index (k)**

  ```cpp
  int ij_to_k(int i, int j, int M);
  ```

  Converts 2D indices on the internal ((M-2)x(M-2)) grid into a 1D index (k).

* **Build the diagonal vectors (a(N)) and (b(N)) from (M), (\Delta t) and (V)**

  ```cpp
  void build_ab_vectors(int M, double dt, const arma::mat& V,
                        arma::cx_vec& a, arma::cx_vec& b);
  ```

  Given the potential matrix (V) and time step (dt), constructs the complex diagonal entries that encode both the kinetic and potential contributions for all internal grid points.

* **Build sparse matrices (A) and (B) for the 2D Crank–Nicolson scheme**

  ```cpp
  void construct_AB_matrices(const arma::cx_vec& a,
                             const arma::cx_vec& b,
                             std::complex<double> r,
                             int M,
                             arma::sp_cx_mat& A,
                             arma::sp_cx_mat& B);
  ```

  Uses the diagonal vectors `a` and `b` and the parameter `r = i Δt / (2 h^2)` to assemble the sparse matrices (A) and (B) that appear in the linear system (A u^{n+1} = B u^n).

* **Print the structure of a sparse matrix to screen**

  ```cpp
  void print_sp_matrix_structure(const arma::sp_cx_mat& A);
  ```

  Prints a textual representation of the non-zero pattern of a sparse matrix. Useful to verify that the stencil structure of (A) or (B) is correct.

* **Initialise the potential (V(x,y)) for a vertical wall with slits**

  ```cpp
  void init_potential(arma::mat& V,
                      int M,
                      double v0,
                      double wall_thickness_x,
                      double wall_x_pos,
                      double wall_sep_length,
                      double slit_aperture,
                      int n_slits);
  ```

  Fills the real-valued potential matrix (V) according to a wall centred at `wall_x_pos` with thickness `wall_thickness_x` and a chosen number of slits, each of aperture `slit_aperture` and separated by `wall_sep_length`.

* **Print a representation of the potential (V) in the terminal**

  ```cpp
  void print_potential_structure(const arma::mat& V, double threshold);
  ```

  Shows a coarse ASCII “map” of the potential, marking barrier regions and free regions. Helpful to quickly verify the slit geometry.

* **Save a grayscale image of the potential (V) to a PGM file**

  ```cpp
  void save_potential_image(const arma::mat& V,
                            const std::string& filename);
  ```

  Writes `V` as a PDF image.

* **Set up the initial state (u(x,y,t=0)) as a normalised Gaussian wave packet**

  ```cpp
  void init_gaussian_packet(arma::cx_mat& U0,
                            int M,
                            double xc, double yc,
                            double sigma_x, double sigma_y,
                            double p_x, double p_y);
  ```

  Builds the complex Gaussian wave packet on the (M\times M) grid, applies Dirichlet boundary conditions, and normalises it such that the sum of the probability is 1.

* **Extract internal values of (U) ((MxM)) into a vector (u)**

  ```cpp
  arma::cx_vec pack_internal_to_vec(const arma::cx_mat& U, int M);
  ```

  Reads the internal grid points from the matrix (U) and packs them into the 1D vector `u` using the `ij_to_k` mapping. Boundary values (Dirichlet zeros) are excluded.

* **Perform the matrix multiplication (B u^n = b)**

  ```cpp
  void matrix_mult(const arma::sp_cx_mat& B,
                   const arma::cx_vec& u_n,
                   arma::cx_vec& b);
  ```

  Computes the right-hand side vector (b = B u^n) before solving the Crank–Nicolson system.

* **Perform one Crank–Nicolson time step**

  ```cpp
  bool cn_step(const arma::sp_cx_mat& A,
               const arma::sp_cx_mat& B,
               const arma::cx_vec& u,
               arma::cx_vec& u_next);
  ```

  Wraps the two-stage update
  (b = B u) and (A u_next = b) using `arma::spsolve`. Returns `true` if the sparse solve succeeds. This is the core routine used inside the main time-stepping loop.

* **Fill an (MxM) matrix (U) from the internal vector (u)**

  ```cpp
  void unpack_vec_to_internal(arma::cx_mat& U,
                              const arma::cx_vec& u,
                              int M);
  ```

  Performs the inverse of `pack_internal_to_vec`: starts from a zero (MxM) matrix, writes the internal entries from `u` using `ij_to_k`, and leaves boundary points at zero. This is used to reconstruct 2D fields from the 1D state vector for plotting.

---

## 7. Extending the project

The structure of the code makes it straightforward to:

* Change the potential geometry (e.g. other barrier shapes, time-dependent potentials).
* Experiment with different initial states.

```
```

