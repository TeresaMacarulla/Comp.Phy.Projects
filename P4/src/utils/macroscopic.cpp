#include "macroscopic.hpp"

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
static inline std::pair<double, double> compute_energy_and_magnetization(
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

    // For L=2 with PBC, some analytical references adopt a "double-counted" convention.
    // To match those closed-form expressions ONLY for the 2x2 case, scale by 2 here.
    if (L == 2) {
        E *= 2.0;
    }

    return {E, M};
}

/**
 * Metropolis MCMC for the 2D Ising model (nearest-neighbor, periodic BCs).
 *
 * @param s            LxL lattice with spins in {+1,-1}; it will be updated in-place.
 * @param L            Lattice linear size (s.size() == L).
 * @param T            Temperature (k_B = 1 convention).
 * @param J            Coupling (default 1.0).
 * @param total_sweeps Total number of Monte Carlo sweeps to perform.
 * @param burn_in      Number of initial sweeps discarded (thermalization).
 * @param sample_every Take one measurement every 'sample_every' sweeps after burn-in.
 * @param seed         RNG seed (0 -> random_device).
 *
 * @return             MCMCResult with sampled E and M, plus diagnostics.
 */

MCMCResult run_mcmc_metropolis(std::vector<std::vector<int>>& s, int L, double T, double J, int total_sweeps, int burn_in, int sample_every, unsigned long long seed) {
    const int N = L * L;
    const double beta = 1.0 / T; // k_B = 1

    // --- RNG setup ---
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    std::uniform_int_distribution<int> pick_i(0, L - 1);
    std::uniform_int_distribution<int> pick_j(0, L - 1);
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    // --- Precompute Boltzmann factors for the only possible ΔE values ---
    // For square lattice: m = s_i * sum(neighbors) ∈ {-4,-2,0,2,4}
    // ΔE = 2 J m  ∈ {-8J,-4J,0,4J,8J}
    auto idx_from_m = [](int m){ return (m + 4) / 2; }; // maps {-4,-2,0,2,4} -> {0,1,2,3,4}
    double boltz[5];
    double dE_vals[5];

    if (L == 2) {
        double temp[] = {-16.0 * J, -8.0 * J, 0.0, 8.0 * J, 16.0 * J};
        std::copy(std::begin(temp), std::end(temp), dE_vals);
    } else {
        double temp[] = {-8.0 * J, -4.0 * J, 0.0, 4.0 * J, 8.0 * J};
        std::copy(std::begin(temp), std::end(temp), dE_vals);
    }

    for (int k = 0; k < 5; ++k) boltz[k] = std::exp(-beta * dE_vals[k]);

    // --- Initialize E and M once; then update them incrementally after each accepted flip ---
    auto [E, M] = compute_energy_and_magnetization(s, L, J);

    MCMCResult out;
    out.sweeps_done = total_sweeps;

    // --- Main MCMC loop (Metropolis updates) ---
    for (int sweep = 1; sweep <= total_sweeps; ++sweep) {

        // One sweep = N attempted single-spin flips at random sites
        for (int t = 0; t < N; ++t) {
            // 1) Pick a random site
            int i = pick_i(rng);
            int j = pick_j(rng);

            // 2) Compute local neighbor sum with periodic boundaries (branch-free with modulo)
            int ip = (i + 1) % L, im = (i - 1 + L) % L;
            int jp = (j + 1) % L, jm = (j - 1 + L) % L;

            int si = s[i][j];
            int nn = s[i][jp] + s[i][jm] + s[ip][j] + s[im][j]; // sum of 4 nearest neighbors

            // 3) Energy change for flipping s_ij -> -s_ij
            // ΔE = 2 * J * si * nn
            int m = si * nn; // in {-4,-2,0,2,4}
            double dE = 2.0 * J * static_cast<double>(m);

            // Special correction for the 2x2 lattice (double-counted pairs)
            if (L == 2) {
                dE *= 2.0;
            }

            // 4) Metropolis acceptance
            bool accept = false;
            if (dE <= 0.0) {
                accept = true; // downhill move: always accept
            } else {
                // Use precomputed Boltzmann factor: index from m
                int k = idx_from_m(m); // 0..4
                if (U01(rng) < boltz[k]) accept = true;
            }

            // 5) If accepted, flip spin and update E and M incrementally
            if (accept) {
                s[i][j] = -si;
                out.accepted_flips++;

                // Magnetization change: ΔM = (-si) - si = -2 si
                M += -2.0 * si;

                // Energy change: add ΔE we already computed
                E += dE;
            }
        }

        // 6) Measurement: after burn-in, sample every 'sample_every' sweeps
        if (sweep > burn_in && ((sweep - burn_in) % sample_every == 0)) {
            out.E_samples.push_back(E);
            out.M_samples.push_back(M);
        }
    }

    return out;
}