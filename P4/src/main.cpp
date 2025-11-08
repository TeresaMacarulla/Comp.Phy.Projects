#include "initialize_lattice.hpp"

#include <iostream>

// COMPILING AND LINKING FROM P4 FILE: g++ -std=c++17 -O2 -I include/*.hpp src/main.cpp src/utils/*.cpp -larmadillo -o src/main.exe

int main() {
    int L = 10;
    double J = 1.0;

    // Create lattice
    auto lattice = initialize_lattice(L);

    // Compute energy and magnetization
    //auto [E, M] = compute_energy_and_magnetization(lattice, L, J);

    //std::cout << "Energy = " << E << "\n";
    //std::cout << "Magnetization = " << M << "\n";

    for (const auto& row : lattice) {
        for (const auto& val : row) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}