#include "PenningTrap.hpp"
#include <armadillo>
#include <fstream>
#include <iomanip>

// Compiling and linking: g++ -std=c++17 -O2 -I include/utils.h src/main.cpp src/utils.cpp/*.cpp -larmadillo -o src/main.exe

// Function declarations
void evolution_rz(PenningTrap& trap, double T_us, int Nstep);

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

    //Add particles to the trap
    trap.add_particle(p1);
    //trap.add_particle(p2);

    // Simulate a single particle and track z-direction movement
    evolution_rz(trap, 50.0, 100);

    std::cout << "Particle p1:" << std::endl;
    std::cout << "  Charge (q): " << p1.q << std::endl;
    std::cout << "  Mass (m): "   << p1.m << std::endl;
    std::cout << "  Position (r): " << p1.r.t();  // .t() → imprime como fila
    std::cout << "  Velocity (v): " << p1.v.t();

    std::cout << "\nPenning Trap:" << std::endl;
    std::cout << "  B0 = " << trap.B0 << std::endl;
    std::cout << "  V0 = " << trap.V0 << std::endl;
    std::cout << "  d  = " << trap.d  << std::endl;
    std::cout << "  Number of particles: " << trap.particles.size() << std::endl;

    for (std::size_t i = 0; i < trap.particles.size(); ++i) {
    std::cout << "\nParticle " << i << ":" << std::endl;
    std::cout << "  q = " << trap.particles[i].q << std::endl;
    std::cout << "  m = " << trap.particles[i].m << std::endl;
    std::cout << "  r = " << trap.particles[i].r.t();
    std::cout << "  v = " << trap.particles[i].v.t();
    }

    return 0;
}

void evolution_rz(PenningTrap& trap, double T_us, int Nstep)
{
    // Compute RK4 time step h = T/N
    const double h = T_us / Nstep;

    // Open output file
    std::ofstream out("../scripts/rz.txt");
    out << std::scientific << std::setprecision(8);

    // Write initial state (t = 0)
    double t = 0.0;
    out << t << "  " << trap.particles[0].r(2) << "\n";

    // Main RK4 loop
    for (int k = 0; k < Nstep; ++k) {
        trap.evolve_forward_Euler(h);        // advance all particles by h
        t += h;                    // update time (in microseconds)
        out << t << "  " << trap.particles[0].r(2) << "\n";  // log r_z of particle 0
    }

    out.close();
}