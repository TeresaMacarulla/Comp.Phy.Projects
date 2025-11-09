#include "macroscopic.hpp"

#include <iostream>
#include <iomanip>

// COMPILING AND LINKING FROM P4 FILE: g++ -std=c++17 -O2 -I include src/main.cpp src/utils/*.cpp -larmadillo -o src/main.exe

int main() {

    int option;

    // Show message in the terminal
    std::cout << "\n=== ISING MODEL SIMULATION MENU ===\n"
          << " Choose an option\n"  
          << " 1. Apply the Markov Chain Monte Carlo approach to sample spin configurations and compute mean energy and magnetisation, heat capacity and susceptibility per spin for an fixed temperature\n"
          << " 2. \n"
          << " 3. \n"
          << " 4. \n"
          << " 5.\n"
          << " 6. \n";

    std::cin >> option;

    if (option==1 ){

        int L;

        std::cout << "Choose L for the lattice dimension (LxL)\n";
        std::cin >> L;

        double T;

        std::cout << "Choose temperature (J/kB units)\n";
        std::cin >> T;

        int MC_sweeps;
        std::cout << "Choose number of Monte Carlo sweeps\n";
        std::cin >> MC_sweeps;

        int J = 1;

        // Create lattice
        auto lattice = initialize_lattice(L);

        auto res = run_mcmc_metropolis(lattice, L, T, J, MC_sweeps);

        // Post-process results from "run_mcmc_metropolis"
        const int N = L*L;
        double E_mean = 0.0, E2_mean = 0.0, Mabs_mean = 0.0, M2_mean = 0.0;
        for (size_t k = 0; k < res.E_samples.size(); ++k) {
            double E = res.E_samples[k];
            double M = res.M_samples[k];
            E_mean  += E;
            E2_mean += E*E;
            Mabs_mean += std::abs(M);
            M2_mean += M*M;
        }
        int S = (int)res.E_samples.size();
        E_mean  /= S;  E2_mean /= S;  Mabs_mean /= S;  M2_mean /= S;

        double beta = 1.0 / T;
        double Cv_per_spin  = beta*beta * (E2_mean - E_mean*E_mean) / N;
        double chi_per_spin = beta * (M2_mean - Mabs_mean*Mabs_mean) / N;
        double eps_mean = E_mean / N;
        double mags_mean = Mabs_mean / N;

        std::cout << std::fixed << std::setprecision(12);

        std::cout << " For a temperature of T = " << T << "J/kB, the results are : " << "\n"
                  << "Cv/N = " << Cv_per_spin << "\n"
                  << "χ/N = " << chi_per_spin << "\n"
                  << "<ε> = " << eps_mean << "\n"
                  << "<|m|> = " << mags_mean << "\n";
    }

    return 0;
}