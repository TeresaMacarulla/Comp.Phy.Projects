#include "functions.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>

using namespace std;

// COMPILING AND LINKING FROM P4 FILE: g++ -std=c++17 -O2 -I include src/main.cpp src/utils/*.cpp -larmadillo -o src/main.exe -fopenmp

int main() {

    int option;

    // Show message in the terminal
    cout << "\n=== ISING MODEL SIMULATION MENU ===\n"
          << " Choose an option\n"  
          << " 1. Apply the Markov Chain Monte Carlo approach to sample spin configurations and compute mean energy and magnetisation, heat capacity and susceptibility per spin for an fixed temperature\n"
          << " 2. Study the burn-in time (or equilibration time)\n"
          << " 3. Run a simulation choosing the adecuate burn-in Monte Carlo cycles \n"
          << " 4. Choose some temperature values and run simulations with and without parallelization (OpenMP). Compare the timing.\n"
          << " 5. Run parallelization from T=2.1 to T=2.4 J/kB and from L=40 to L=100\n";

    cin >> option;

    int J = 1;

    if (option == 1){

        int L;
        cout << "Choose L for the lattice dimension (LxL)\n";
        cin >> L;

        double T;
        cout << "Choose temperature (J/kB units)\n";
        cin >> T;

        int MC_sweeps;
        cout << "Choose number of Monte Carlo sweeps\n";
        cin >> MC_sweeps;

        // Create lattice
        auto lattice = initialize_lattice(L);

        auto res = run_mcmc_metropolis(lattice, L, T, J, MC_sweeps);

        auto parameters = compute_parameters (res.E_samples, res.M_samples, L, T);

        cout << fixed << setprecision(12);

        cout << " For a temperature of T = " << T << "J/kB, the results are : " << "\n"
                  << "Cv/N = " << parameters.Cv_per_spin << "\n"
                  << "χ/N = " << parameters.chi_per_spin << "\n"
                  << "<ε> = " << parameters.eps_mean << "\n"
                  << "<|m|> = " << parameters.mags_mean << "\n"
                  << "Accepted flips = " << res.accepted_flips << "\n";
    }

    if (option == 2){

        int L;
        cout << "Choose L for the lattice dimension (LxL)\n";
        cin >> L;

        double T;
        cout << "Choose temperature (J/kB units)\n";
        cin >> T;

        int MC_sweeps;
        cout << "Choose number of Monte Carlo sweeps\n";
        cin >> MC_sweeps;

        // Output file 
        string filename = "txt/epsilon_T_" + to_string(T) + ".txt";
        ofstream out(filename);
        out << scientific << setprecision(12);

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

        cout << "SEE THE RESULTS: run <<equilibration_results.py>> in <<scripts>>\n"
                  << "This script is particularly made for T=1 and T=2.4\n";
    }

    if (option == 3){

        int L;
        cout << "Choose L for the lattice dimension (LxL)\n";
        cin >> L;

        double T;
        cout << "Choose temperature (J/kB units)\n";
        cin >> T;

        int MC_sweeps;
        cout << "Choose number of Monte Carlo sweeps\n";
        cin >> MC_sweeps;

        cout << "Adecuate Burn-in sweeps for termalization:\n"
                  << "T = 1 J/kB --> 100 burn-in\n"
                  << "T = 2.4 J/kB --> 10^5 burn-in\n";
        int burn_in;
        cout << "Choose an adecuate number of Burn-in sweeps for termalization\n";
        cin >> burn_in;

        // Output file 
        string filename = "txt/eps_T_" + to_string(T) + "_burn-in_" + to_string(burn_in) + ".txt";
        ofstream out(filename);
        out << scientific << setprecision(12);

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

        cout << "You should run this option until you have two 'burn-in' .txt files (T=1 and T=2.4 J/kB)\n"
                  << "Then, you can run <<histogram.py>> in P4/scripts to obtain the histograms of the probability of <ε>\n";

    }

    if (option == 4){

        vector<double> temperatures;
        cout << "Enter temperature values in J/kB units (separated by spaces):\n";
        string line;
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        getline(cin, line);      // read whole line
        istringstream iss(line);
        double T;
        while (iss >> T) {
           temperatures.push_back(T);
        }

        int L;
        cout << "Choose L for the lattice dimension (LxL)\n";
        cin >> L;

        int MC_sweeps;
        cout << "Choose number of Monte Carlo sweeps\n";
        cin >> MC_sweeps;

        cout << "Adecuate Burn-in sweeps for termalization:\n"
                  << "T = 1 J/kB --> 100 burn-in\n"
                  << "T = 2.4 J/kB --> 10^5 burn-in\n";
        int burn_in;
        cout << "Choose an adecuate number of Burn-in sweeps for termalization (considering your temperatures)\n";
        cin >> burn_in;

        // Start measuring time (parallel code)
        auto tp1 = chrono::high_resolution_clock::now();
      
        #pragma omp parallel for
        for (double temp : temperatures) {
            auto lat = initialize_lattice(L);
            auto res = run_mcmc_metropolis(lat, L, temp, J, MC_sweeps, burn_in);
        }

        // Stop measuring time (parallel code)
        auto tp2 = chrono::high_resolution_clock::now();

        // Calculate the elapsed time (parallel code)
        double tp = chrono::duration<double>(tp2 - tp1).count();
        cout << "Timing for parallel code: " << tp << "\n";

        // Start measuring time (serial code)
        auto ts1 = chrono::high_resolution_clock::now();

        for (double temp : temperatures) {
            auto lat = initialize_lattice(L);
            auto res = run_mcmc_metropolis(lat, L, temp, J, MC_sweeps, burn_in);
        }

        // Stop measuring time (serial code)
        auto ts2 = chrono::high_resolution_clock::now();

        // Calculate the elapsed time (serial code)
        double ts = chrono::duration<double>(ts2 - ts1).count();
        cout << "Timing for serial code: " << ts << "\n";

        double speedup = ts/tp;
        cout << "The speed-up factor is: " << speedup << "\n";
    }

    if (option == 5){

        vector<int> L_size = {40, 60, 80, 100};
      
        for (int L : L_size) {
            const int N = L*L;

            // Output file 
            string filename = "txt/parameters_L_" + to_string(L) + ".txt";
            ofstream out(filename);
            out << scientific << setprecision(12);

            vector<double> temperatures;
            // 2.10 → 2.19
            for (double T = 2.10; T < 2.20; T += 0.007)
                temperatures.push_back(T);
            // 2.20 → 2.35
            for (double T = 2.20; T < 2.35 + 1e-9; T += 0.0043)
                temperatures.push_back(T);
            // 2.36 → 2.40
            for (double T = 2.36; T < 2.41; T += 0.007)
                temperatures.push_back(T);

            #pragma omp parallel for
            for (int i = 0; i < temperatures.size(); ++i) {
                double temp = temperatures[i];
                auto lat = initialize_lattice(L);
                auto res = run_mcmc_metropolis(lat, L, temp, J, 100000, 1000);

                auto parameters = compute_parameters (res.E_samples, res.M_samples, L, temp);
                out << temp << " " << parameters.Cv_per_spin << " " << parameters.chi_per_spin << " " << parameters.eps_mean << " " << parameters.mags_mean << "\n";
            }
            out.close();
        }

        cout << "Find the results in P4/txt fonder, in the files 'parameters_L_40' and so on." << "\n";

    }

    return 0;

}