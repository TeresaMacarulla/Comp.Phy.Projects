#include "PenningTrap.hpp"
#include "PenningTrap_Vt.hpp"
#include "Evolution.hpp"

// COMPILING AND LINKING : g++ -std=c++17 -O2 -I include/utils.h src/main.cpp src/utils.cpp/*.cpp -larmadillo -o src/main.exe

// std::cout << " rz " << r(2) << " Total Force " << F(2) << "\n";

int main() {

    // Define B0, V0 and d in the trap
    PenningTrap trap(96.5, 2.41e6, 500.0, {});

    // Define B0, V0 and d in the trap.
    PenningTrap_Vt trap_Vt(96.5, 2.41e6, 500.0, 0.1, 0.2, {});

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

    // Make a copy of trap
    //PenningTrap trap1 = trap;

    // Simulate a single particle and track z-direction movement
    //evolution_rz(trap1, 50.0, 1000);

    // Simulate two particles and track rx and ry movement 
    //trap1 = trap;
    //evolution_two(trap1, 50.0, 1000, 0); // last input =0 (no interaction) or =1 (interaction)

    // Compare analytical solution, RK4 and Euler for single particle
    //trap1 = trap;
    //evolution_nsteps(4000, trap1);

    // Fill trap_Vt with particles with randomly generated initial positions and velocities
    fill_trap(trap_Vt, 0.1, 0.2, 100);
  
    // Try f amplitudes {0.1, 0.4, 0.7} for V0(y)
    std::vector<double> f = {0.1, 0.4, 0.7};

    // Loop through f elements
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
            std::size_t n_inside = evolution_many(trap_Vt, 500.0, 10000, 0);

            out << omegav << " " << n_inside << " " << "\n";

        }

        out.close();

    }

    return 0;
}

