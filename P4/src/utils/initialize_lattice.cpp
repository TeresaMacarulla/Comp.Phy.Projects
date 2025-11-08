#include "initialize_lattice.hpp"

// Function 1: Create a lattice LxL with random spin values (+1 or -1)
std::vector<std::vector<int>> initialize_lattice(int L, double p_up, unsigned int seed) {
    std::vector<std::vector<int>> lattice(L, std::vector<int>(L));
    std::mt19937 rng(seed ? seed : std::random_device{}());
    std::bernoulli_distribution spin_up(p_up);

    // Assign +1 or -1 to each lattice site
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            lattice[i][j] = spin_up(rng) ? +1 : -1;
        }
    }
    return lattice;
}

// Function 2: Compute total Energy (E) and Magnetization (M) of the lattice
std::pair<double, double> compute_energy_and_magnetization(
        const std::vector<std::vector<int>>& lattice, int L, double J) {
    
    double E = 0.0;
    double M = 0.0;

    // Loop over all lattice sites
    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            int s = lattice[i][j];
            M += s; // accumulate total magnetization

            // Periodic boundary conditions: each spin interacts with right and down neighbors
            int right = lattice[i][(j + 1) % L];
            int down  = lattice[(i + 1) % L][j];

            // Add nearest-neighbor contributions
            E += -J * s * (right + down);
        }
    }

    // E is the total energy; M is the total magnetization
    return {E, M};
}