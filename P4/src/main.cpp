#include "functions.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>

// COMPILING AND LINKING FROM P4 FILE: g++ -std=c++17 -O2 -I include src/main.cpp src/utils/*.cpp -larmadillo -o src/main.exe -fopenmp

int main() {

    int option;

    // Show message in the terminal
    std::cout << "\n=== ISING MODEL SIMULATION MENU ===\n"
          << " Choose an option\n"  
          << " 1. Apply the Markov Chain Monte Carlo approach to sample spin configurations and compute mean energy and magnetisation, heat capacity and susceptibility per spin for an fixed temperature\n"
          << " 2. Study the burn-in time (or equilibration time)\n"
          << " 3. Run a simulation choosing the adecuate burn-in Monte Carlo cycles \n"
          << " 4. Choose some temperature values and run simulations with and without parallelization (OpenMP). Compare the timing.\n"
          << " 5. Run parallelization from T=2.1 to T=2.4 J/kB and from L=40 to L=100\n"
          << " 6. \n";

    std::cin >> option;

    int J = 1;

    if (option == 1){

        int L;
        std::cout << "Choose L for the lattice dimension (LxL)\n";
        std::cin >> L;

        double T;
        std::cout << "Choose temperature (J/kB units)\n";
        std::cin >> T;

        int MC_sweeps;
        std::cout << "Choose number of Monte Carlo sweeps\n";
        std::cin >> MC_sweeps;

        // Create lattice
        auto lattice = initialize_lattice(L);

        auto res = run_mcmc_metropolis(lattice, L, T, J, MC_sweeps);

        auto parameters = compute_parameters (res.E_samples, res.M_samples, L, T);

        std::cout << std::fixed << std::setprecision(12);

        std::cout << " For a temperature of T = " << T << "J/kB, the results are : " << "\n"
                  << "Cv/N = " << parameters.Cv_per_spin << "\n"
                  << "χ/N = " << parameters.chi_per_spin << "\n"
                  << "<ε> = " << parameters.eps_mean << "\n"
                  << "<|m|> = " << parameters.mags_mean << "\n"
                  << "Accepted flips = " << res.accepted_flips << "\n";
    }

    if (option == 2){

        int L;
        std::cout << "Choose L for the lattice dimension (LxL)\n";
        std::cin >> L;

        double T;
        std::cout << "Choose temperature (J/kB units)\n";
        std::cin >> T;

        int MC_sweeps;
        std::cout << "Choose number of Monte Carlo sweeps\n";
        std::cin >> MC_sweeps;

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

        std::cout << "SEE THE RESULTS: run <<equilibration_results.py>> in <<scripts>>\n"
                  << "This script is particularly made for T=1 and T=2.4\n";
    }

    if (option == 3){

        int L;
        std::cout << "Choose L for the lattice dimension (LxL)\n";
        std::cin >> L;

        double T;
        std::cout << "Choose temperature (J/kB units)\n";
        std::cin >> T;

        int MC_sweeps;
        std::cout << "Choose number of Monte Carlo sweeps\n";
        std::cin >> MC_sweeps;

        std::cout << "Adecuate Burn-in sweeps for termalization:\n"
                  << "T = 1 J/kB --> 100 burn-in\n"
                  << "T = 2.4 J/kB --> 10^5 burn-in\n";
        int burn_in;
        std::cout << "Choose an adecuate number of Burn-in sweeps for termalization\n";
        std::cin >> burn_in;

        // Output file 
        std::string filename = "txt/eps_T_" + std::to_string(T) + "_burn-in_" + std::to_string(burn_in) + ".txt";
        std::ofstream out(filename);
        out << std::scientific << std::setprecision(12);

        // Create lattices and run Monte Carlo cycles
        auto lat = initialize_lattice(L);
        auto res = run_mcmc_metropolis(lat, L, T, J, MC_sweeps, burn_in);

        // Post-process results from "run_mcmc_metropolis"
        const int N = L*L;
        double eps = 0.0;

        //POSSIBLE PARALLELIZING
        for (size_t k = 0; k < res.E_samples.size(); ++k) {
            double E = res.E_samples[k];
            eps = E/N;
 
            out << (k+1) << " " << eps << " " << E << " " << "\n";
        }
        out.close();

        std::cout << "You should run this option until you have two 'burn-in' .txt files (T=1 and T=2.4 J/kB)\n"
                  << "Then, you can run <<histogram.py>> in P4/scripts to obtain the histograms of the probability of <ε>\n";

    }

    if (option == 4){

        std::vector<double> temperatures;
        std::cout << "Enter temperature values in J/kB units (separated by spaces):\n";
        std::string line;
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::getline(std::cin, line);      // read whole line
        std::istringstream iss(line);
        double T;
        while (iss >> T) {
           temperatures.push_back(T);
        }

        int L;
        std::cout << "Choose L for the lattice dimension (LxL)\n";
        std::cin >> L;

        int MC_sweeps;
        std::cout << "Choose number of Monte Carlo sweeps\n";
        std::cin >> MC_sweeps;

        std::cout << "Adecuate Burn-in sweeps for termalization:\n"
                  << "T = 1 J/kB --> 100 burn-in\n"
                  << "T = 2.4 J/kB --> 10^5 burn-in\n";
        int burn_in;
        std::cout << "Choose an adecuate number of Burn-in sweeps for termalization (considering your temperatures)\n";
        std::cin >> burn_in;

        // Start measuring time (parallel code)
        auto tp1 = std::chrono::high_resolution_clock::now();
      
        #pragma omp parallel for
        for (double temp : temperatures) {
            auto lat = initialize_lattice(L);
            auto res = run_mcmc_metropolis(lat, L, temp, J, MC_sweeps, burn_in);
        }

        // Stop measuring time (parallel code)
        auto tp2 = std::chrono::high_resolution_clock::now();

        // Calculate the elapsed time (parallel code)
        double tp = std::chrono::duration<double>(tp2 - tp1).count();
        std::cout << "Timing for parallel code: " << tp << "\n";


        // Start measuring time (serial code)
        auto ts1 = std::chrono::high_resolution_clock::now();

        for (double temp : temperatures) {
            auto lat = initialize_lattice(L);
            auto res = run_mcmc_metropolis(lat, L, temp, J, MC_sweeps, burn_in);
        }

        // Stop measuring time (serial code)
        auto ts2 = std::chrono::high_resolution_clock::now();

        // Calculate the elapsed time (serial code)
        double ts = std::chrono::duration<double>(ts2 - ts1).count();
        std::cout << "Timing for serial code: " << ts << "\n";

        double speedup = ts/tp;
        std::cout << "The speed-up factor is: " << speedup << "\n";

    }

    if (option == 5){

        std::vector<int> L_size = {40, 60, 80, 100};
      
        for (int L : L_size) {
            const int N = L*L;

            // Output file 
            std::string filename = "txt/parameters_L_" + std::to_string(L) + ".txt";
            std::ofstream out(filename);
            out << std::scientific << std::setprecision(12);

            #pragma omp parallel for
            for (int T = 210; T < 241; T += 1) {
                double temp = T/100.0;
                auto lat = initialize_lattice(L);
                auto res = run_mcmc_metropolis(lat, L, temp, J, 100000, 1000);

                auto parameters = compute_parameters (res.E_samples, res.M_samples, L, T);
                out << temp << " " << parameters.Cv_per_spin << " " << parameters.chi_per_spin << " " << parameters.eps_mean << " " << parameters.mags_mean << "\n";

                /*
                std::cout << std::fixed << std::setprecision(12);

                std::cout << " For a temperature of T = " << T << "J/kB, the results are : " << "\n"
                  << "Cv/N = " << parameters.Cv_per_spin << "\n"
                  << "χ/N = " << parameters.chi_per_spin << "\n"
                  << "<ε> = " << parameters.eps_mean << "\n"
                  << "<|m|> = " << parameters.mags_mean << "\n"
                  << "Accepted flips = " << res.accepted_flips << "\n";
                */
            }
            out.close();
        }

    }

    if (option == 6){

    }

    return 0;
}