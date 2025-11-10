#include "macroscopic.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>

// COMPILING AND LINKING FROM P4 FILE: g++ -std=c++17 -O2 -I include src/main.cpp src/utils/*.cpp -larmadillo -o src/main.exe

int main() {

    int option;

    // Show message in the terminal
    std::cout << "\n=== ISING MODEL SIMULATION MENU ===\n"
          << " Choose an option\n"  
          << " 1. Apply the Markov Chain Monte Carlo approach to sample spin configurations and compute mean energy and magnetisation, heat capacity and susceptibility per spin for an fixed temperature\n"
          << " 2. Study the burn-in time (or equilibration time)\n"
          << " 3. \n"
          << " 4. \n"
          << " 5.\n"
          << " 6. \n";

    std::cin >> option;

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

    if (option==1 ){

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
                  << "<|m|> = " << mags_mean << "\n"
                  << "Accepted flips = " << res.accepted_flips << "\n";
    }

    if (option==2 ){

        // Output file 
        std::string filename = "txt/epsilon_T_" + std::to_string(T) + ".txt";
        std::ofstream out(filename);
        out << std::scientific << std::setprecision(12);

        // Create lattices (ordered and disordered)
        auto lat_ordered = initialize_lattice(L, true);
        auto lat_disordered = initialize_lattice(L);

        auto res_ordered = run_mcmc_metropolis(lat_ordered, L, T, J, MC_sweeps);
        auto res_disordered = run_mcmc_metropolis(lat_disordered, L, T, J, MC_sweeps);

        // Post-process results from "run_mcmc_metropolis"
        const int N = L*L;
        double eps_o = 0.0;
        double eps_mean_o = 0.0;
        double eps_sum_o = 0.0;

        double eps_d = 0.0;
        double eps_mean_d = 0.0;
        double eps_sum_d = 0.0;


        for (size_t k = 0; k < res_ordered.E.size(); ++k) {
            double E_o = res_ordered.E[k];
            eps_o = E_o/N;
            eps_sum_o += eps_o;
            eps_mean_o = eps_sum_o/(k+1);

            double E_d = res_disordered.E[k];
            eps_d = E_d/N;
            eps_sum_d += eps_d;
            eps_mean_d = eps_sum_d/(k+1);

            out << (k+1) << " " << eps_o << " " << eps_mean_o << " " << eps_d << " " << eps_mean_d << " "<< "\n";
        }
        out.close();

        std::cout << "SEE THE RESULTS: run <<equilibration_results.py>> in <<scripts>>\n";
    }

    return 0;
}