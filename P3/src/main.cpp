#include "PenningTrap.hpp"
#include "PenningTrap_Vt.hpp"
#include "Evolution.hpp"

// COMPILING AND LINKING : g++ -std=c++17 -O2 -I include/utils.h src/main.cpp src/utils.cpp/*.cpp -larmadillo -o src/main.exe

// Not recommendable to run all the code at once, comment the parts that are no needed

int main() {

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

    // Simulate a single particle and track z-direction movement
    PenningTrap trap1 = trap; // Make a copy of trap  
    evolution_rz(trap1, 50.0, 1000);

    // Simulate two particles and track rx and ry movement 
    trap1 = trap;
    evolution_two(trap1, 50.0, 1000, 0); // last input =0 (no interaction) or =1 (interaction)

    // Compare analytical solution, RK4 and Euler for single particle
    trap1 = trap;
    evolution_nsteps(4000, trap1);  // Applies the cases nsteps={4000, 8000, 16000, 32000} all at once


    // V(t) TIME DEPENDENCY IN V0-----------------------------------------------------------------------

    // Define B0, V0 and d in the trap with V(t)
    PenningTrap_Vt trap_Vt(96.5, 2.41e6, 500.0, 0.1, 0.2, {});

    // Fill trap_Vt with particles with randomly generated initial positions and velocities
    fill_trap(trap_Vt, 0.1, 0.2, 100);
  
    // Try f amplitudes {0.1, 0.4, 0.7} for V0(y)
    std::vector<double> f = {0.1, 0.4, 0.7};

    // Loop through f elements: 
    // first frequency range takes ~8 minutes, the zoom part below takes ~5 minutes with no interactions and too much minutes with interactions
    for (std::size_t i = 0; i < f.size(); ++i) {

        trap_Vt.f = f[i];

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

        // Zoom in frequency
        // Output file 
        std::string filename_zoom = "../scripts/trapped_vs_freq_zoom" + std::to_string(f[i]) + ".txt"; 
        //std::string filename_zoom = "../scripts/trap_vs_freq_zoom_inter" + std::to_string(f[i]) + ".txt"; 
        std::ofstream out_zoom(filename_zoom);
        out_zoom << std::scientific << std::setprecision(12);

        for (double omegav = 1.0; omegav <= 1.75; omegav += 0.01) {
    
            trap_Vt.omegav = omegav;
            PenningTrap_Vt trap_Vt1 = trap_Vt;

            // Simulate 100 particles during 500 microseconds
            evolution_many(trap_Vt1, 500.0, 10000, 0);  // change the last input: =0 for no interactions =1 for interactions. 
            std::size_t n_inside = trap_Vt1.inside_particles();

            out_zoom << omegav << " " << n_inside << " " << "\n";

        }

        out_zoom.close();

    }

    return 0;
}

