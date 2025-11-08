#ifndef __INITIALIZE_LATTICE_HPP__
#define __INITIALIZE_LATTICE_HPP__

#include <vector>
#include <random>
#include <utility>   // for std::pair

// Function 1: Create a lattice LxL with random spin values (+1 or -1)
std::vector<std::vector<int>> initialize_lattice(int L, double p_up = 0.5, unsigned int seed = 0);

// Function 2: Compute total Energy (E) and Magnetization (M) of the lattice
std::pair<double, double> compute_energy_and_magnetization(const std::vector<std::vector<int>>& lattice, int L, double J = 1.0);

#endif