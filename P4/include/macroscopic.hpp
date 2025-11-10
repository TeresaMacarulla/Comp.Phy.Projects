#ifndef __MACROSCOPIC_HPP__
#define __MACROSCOPIC_HPP__

#include <vector>
#include <random>
#include <utility>   
#include <cmath>

// Function 1: Create a lattice LxL with random spin values (+1 or -1)
std::vector<std::vector<int>> initialize_lattice(int L, double p_up = 0.5, unsigned int seed = 0, bool ordered = false);

// Function 2: Compute total Energy (E) and Magnetization (M) of the lattice
static inline std::pair<double, double> compute_energy_and_magnetization(const std::vector<std::vector<int>>& lattice, int L, double J = 1.0);

// Function 3: MCMC rutine

struct MCMCResult {std::vector<double> E; std::vector<double> E_samples; std::vector<double> M_samples; int sweeps_done = 0; long long accepted_flips = 0;};

MCMCResult run_mcmc_metropolis(std::vector<std::vector<int>>& s, int L, double T, double J,
    int total_sweeps,
    int burn_in = 500, //1000,
    int sample_every = 60,
    unsigned long long seed = 0ULL);


#endif