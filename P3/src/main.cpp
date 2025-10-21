#include "PenningTrap.hpp"
#include "PenningTrap_Vt.hpp"
#include "Evolution.hpp"
#include <string>
#include <iostream>

// COMPILING AND LINKING FROM P3 FILE: g++ -std=c++17 -O2 -I include/utils.h src/main.cpp src/utils.cpp/*.cpp -larmadillo -o src/main.exe
// RUNNING FROM P3 FILE: ./src/main.exe
// Not recommendable to run all the code at once, comment the parts that are no needed

int main() {

    int option;

    // Show message in the terminal
    std::cout << "\n=== PENNING TRAP SIMULATION MENU ===\n"
          << " 1. Single particle: track z-motion\n"
          << " 2. Two particles: track x–y motion\n"
          << " 3. Compare analytical, RK4, and Euler methods\n"
          << " 4. V0(t): scan f={0.1,0.4,0.7}, ω∈[0.2,2.5], no interactions "
             "(This takes ~3 minutes)\n"
          << " 5. V0(t): zoom f={0.1,0.4,0.7}, ω∈[1.0,1.74], no interactions "
             "(This takes ~2 minutes)\n"
          << " 6. V0(t): zoom f={0.1,0.4,0.7}, ω∈[1.0,1.74], with interactions "
             "(WARNING: very long runtime)\n\n";

    std::cin >> option;

    // **FIRST PART: V0 TIME INDEPENDENT** -----------------------------------------------------------------------

    if (option==1 || option==2 || option==3){
        
        // Define B0, V0 and d in the trap
        PenningTrap trap(96.5, 2.41e6, 500.0, {});

        // Define charge, mass, initial r and initial v of both particles (Ca+)
        arma::vec r1 = {20.0, 0.0, 20.0};
        arma::vec v1 = {0.0, 25.0, 0.0};
        Particle p1(1.0, 40.1, r1, v1);
        arma::vec r2 = {25.0, 25.0, 0.0};
        arma::vec v2 = {0.0, 40.0, 5.0};
        Particle p2(1.0, 40.1, r2, v2);

        // Add particles to the traps
        trap.add_particle(p1);
        trap.add_particle(p2);

        if (option==1){

            // Simulate a single particle and track z-direction movement
            PenningTrap trap1 = trap; // Make a copy of trap  
            evolution_rz(trap1, 50.0, 1000);

            std::cout << "Data stored in /scripts/rz.txt. Run /scripts/rz_vs_t.py python file" << "\n" <<
            " Do you want to perform another action? Choose between [2, 6]" << "\n" <<
            " Do you want to exit? Choose 0" << "\n";
            std::cin >> option;       

            if (option==0){
                std::cout << "Exiting program...\n";
                return 0;
            }

        }

        if (option==2){

            // Simulate two particles and track rx and ry movement 
            PenningTrap trap1 = trap; // Make a copy of trap 
            evolution_two(trap1, 50.0, 1000, 0); // last input =0 (no interaction) or =1 (interaction)

            std::cout << "Data stored in /scripts/ry_rx_nointer.txt, ry_rx_yesinter.txt, rx_vx_nointer.txt, etc. Run /scripts/two_particles.py python file" << "\n" <<
            " Do you want to perform another action? Choose between [3, 6]" << "\n" <<
            " Do you want to exit? Choose 0" << "\n";
            std::cin >> option; 
        
            if (option==0){
                std::cout << "Exiting program...\n";
                return 0;
            }

        }


        if (option==3){

            // Compare analytical solution, RK4 and Euler for single particle
            PenningTrap trap1 = trap; // Make a copy of trap 
            evolution_nsteps(4000, trap1);  // Applies the cases nsteps={4000, 8000, 16000, 32000} all at once

            std::cout << "Data stored in /scripts/r_error_4000.txt, r_error_8000.txt, r_error_16000.txt, r_error_32000.txt. Run /scripts/r_error.py python file" << "\n" <<
            " Do you want to perform another action? Choose between [4, 6]" << "\n" <<
            " Do you want to exit? Choose 0" << "\n";
            std::cin >> option; 
        
            if (option==0){
                std::cout << "Exiting program...\n";
                return 0;
            }
        }
    }

    // **SECOND PART: V(t) TIME DEPENDENCY IN V0**-----------------------------------------------------------------------

    if (option==4 || option==5 || option==6){

        // Define B0, V0 and d in the trap with V(t)
        PenningTrap_Vt trap_Vt(96.5, 2.41e6, 500.0, 0.1, 0.2, {});

        // Fill trap_Vt with particles with randomly generated initial positions and velocities
        fill_trap(trap_Vt, 0.1, 0.2, 100);
  
        // Try f amplitudes {0.1, 0.4, 0.7} for V0(y)
        std::vector<double> f = {0.1, 0.4, 0.7};

        // Loop through f elements: 
        // first frequency range takes ~8 minutes, the zoom part below takes ~5 minutes with no interactions and too much time with interactions
    
        for (std::size_t i = 0; i < f.size(); ++i) {

            trap_Vt.f = f[i];

            if (option==4){
                // Output file 
                std::string filename = "../scripts/trapped_vs_frequency" + std::to_string(f[i]) + ".txt";
                std::ofstream out(filename);
                out << std::scientific << std::setprecision(12);

                for (double omegav = 0.2; omegav <= 2.5; omegav += 0.02) {
    
                    trap_Vt.omegav = omegav;
                    PenningTrap_Vt trap_Vt1 = trap_Vt;

                    // Simulate 100 particles during 500 microseconds
                    evolution_many(trap_Vt1, 500.0, 10000, 0);
                    std::size_t n_inside = trap_Vt1.inside_particles();

                    out << omegav << " " << n_inside << " " << "\n";
                }
                out.close();

                std::cout << "Data stored in /scripts/trapped_vs_frequency0.100000, trapped_vs_frequency0.400000, trapped_vs_frequency0.700000. Run /scripts/trapped_vs_frequency.py python file" << "\n" <<
                " Do you want to perform another action? Choose between [5, 6]" << "\n" <<
                " Do you want to exit? Choose 0" << "\n";
                std::cin >> option;

                if (option==0){
                    std::cout << "Exiting program...\n";
                    return 0;
                }
            }

            if (option==5){            
                // Zoom in frequency
                // Output file 
                std::string filename_zoom = "../scripts/trapped_vs_freq_zoom" + std::to_string(f[i]) + ".txt";  
                std::ofstream out_zoom(filename_zoom);
                out_zoom << std::scientific << std::setprecision(12);

                for (double omegav = 1.0; omegav <= 1.75; omegav += 0.01) {
    
                    trap_Vt.omegav = omegav;
                    PenningTrap_Vt trap_Vt1 = trap_Vt;

                    // Simulate 100 particles during 500 microseconds
                    evolution_many(trap_Vt1, 500.0, 10000, 0);  // last input: =0 no interactions  
                    std::size_t n_inside = trap_Vt1.inside_particles();

                    out_zoom << omegav << " " << n_inside << " " << "\n";
                }
                out_zoom.close();

                std::cout << "Data stored in /scripts/trapped_vs_freq_zoom0.100000, trapped_vs_freq_zoom0.400000, trapped_vs_freq_zoom0.700000. Run /scripts/trapped_vs_frequency.py python file" << "\n" <<
                " Do you want to perform another action? choosing 6 takes to much time to run" << "\n" <<
                " Do you want to exit? Choose 0" << "\n";
                std::cin >> option;

                if (option==0){
                    std::cout << "Exiting program...\n";
                    return 0;
                }   
            }            

            if (option==6){ 
                // Zoom in frequency with interactions
                // Output file 
                std::string filename_zoom = "../scripts/trap_vs_freq_zoom_inter" + std::to_string(f[i]) + ".txt";  
                std::ofstream out_zoom(filename_zoom);
                out_zoom << std::scientific << std::setprecision(12);

                for (double omegav = 1.0; omegav <= 1.75; omegav += 0.01) {
    
                    trap_Vt.omegav = omegav;
                    PenningTrap_Vt trap_Vt1 = trap_Vt;

                    // Simulate 100 particles during 500 microseconds
                    evolution_many(trap_Vt1, 500.0, 10000, 1);  // last input: =1 interactions. 
                    std::size_t n_inside = trap_Vt1.inside_particles();

                    out_zoom << omegav << " " << n_inside << " " << "\n";
                }
                out_zoom.close();

                std::cout << "Data stored in /scripts/trap_vs_freq_zoom_inter0.100000, trap_vs_freq_zoom_inter0.400000, trap_vs_freq_zoom_inter0.700000. Run /scripts/trapped_vs_frequency.py python file" << "\n" <<
                " Do you want to perform another action? Choose between [1, 6]" << "\n" <<
                " Do you want to exit? Choose 0" << "\n";
                std::cin >> option;

                if (option==0){
                    std::cout << "Exiting program...\n";
                    return 0;
                }                   

            }

        }

    }

    return 0;
}

