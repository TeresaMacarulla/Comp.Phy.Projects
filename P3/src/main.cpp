#include "PenningTrap.hpp"
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>

// COMPILING AND LINKING : g++ -std=c++17 -O2 -I include/utils.h src/main.cpp src/utils.cpp/*.cpp -larmadillo -o src/main.exe

// Function declarations
void evolution_rz(PenningTrap& trap, double T_us, int Nstep);
void evolution_two(PenningTrap& trap, double T_us, int Nstep, int inter = 1);
void evolution_RK4_Euler(PenningTrap& trap, double T_us, int Nstep);

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
    trap.add_particle(p2);

    // Simulate a single particle and track z-direction movement
    //evolution_rz(trap, 50.0, 1000);

    // Simulate two particles and track rx and ry movement 
    //evolution_two(trap, 50.0, 1000, 0); // last input =0 (no interaction) or =1 (interaction)

    // Compare analytical solution, RK4 and Euler for single particle
    int i = 4000;
    for (int k = 0; k < 4; ++k) {
        evolution_RK4_Euler(trap, 50.0, i);
        i = i * 2;
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
        trap.evolve_RK4(h); 
        t += h;                    
        out << t << "  " << trap.particles[0].r(2) << "\n";  // log r_z of particle 0
    }

    out.close();
}

void evolution_two(PenningTrap& trap, double T_us, int Nstep, int inter)
{
    // inter = 1 there are interactions between particles
    // inter = 0 there are no interactions between particles

    if (inter != 1 && inter != 0) {
        std::cout << "Error: 'interaction' must be 1 or 0." << std::endl;
        std::exit(EXIT_FAILURE);  
    }

    // Compute RK4 time step h = T/N
    const double h = T_us / Nstep;

    std::ofstream out1;
    std::ofstream out2;
    std::ofstream out3;

    // Open output file depending on 'inter'
    if (inter == 0) {
        out1.open("../scripts/ry_rx_nointer.txt"); // File with ry and rx for both particles (no interactions)
        out2.open("../scripts/rx_vx_nointer.txt"); // File with rx and vx for both particles (no interactions)
        out3.open("../scripts/rz_vz_nointer.txt"); // File with rz and vz for both particles (no interactions)
    } else {
        out1.open("../scripts/ry_rx_yesinter.txt"); // File with ry and rx for both particles (with interactions)
        out2.open("../scripts/rx_vx_yesinter.txt"); // File with rx and vx for both particles (with interactions)
        out3.open("../scripts/rz_vz_yesinter.txt"); // File with rz and vz for both particles (with interactions)
    }

    out1 << std::scientific << std::setprecision(12);
    out2 << std::scientific << std::setprecision(12);
    out3 << std::scientific << std::setprecision(12);

    // Write initial state (t = 0)
    double t = 0.0;
    out1 << t << "  " << trap.particles[0].r(0) << " "<< trap.particles[0].r(1) << " "<< trap.particles[1].r(0) << " "<<trap.particles[1].r(1) << "\n";
    out2 << t << "  " << trap.particles[0].r(0) << " "<< trap.particles[0].v(0) << " "<< trap.particles[1].r(0) << " "<<trap.particles[1].v(0) << "\n";
    out3 << t << "  " << trap.particles[0].r(2) << " "<< trap.particles[0].v(2) << " "<< trap.particles[1].r(2) << " "<<trap.particles[1].v(2) << "\n";

    // Main RK4 loop
    for (int k = 0; k < Nstep; ++k) {
        trap.evolve_RK4(h, inter); 
        t += h;                    
        out1 << t << "  " << trap.particles[0].r(0) << " "<< trap.particles[0].r(1) << " "<< trap.particles[1].r(0) << " "<<trap.particles[1].r(1) << "\n";
        out2 << t << "  " << trap.particles[0].r(0) << " "<< trap.particles[0].v(0) << " "<< trap.particles[1].r(0) << " "<<trap.particles[1].v(0) << "\n";
        out3 << t << "  " << trap.particles[0].r(2) << " "<< trap.particles[0].v(2) << " "<< trap.particles[1].r(2) << " "<<trap.particles[1].v(2) << "\n";  
    }

    out1.close();
    out2.close();
    out3.close();
}

void evolution_RK4_Euler(PenningTrap& trap, double T_us, int Nstep)
{
    // For RK4 evolution
    PenningTrap trap_RK4 = trap;

    // For Euler evolution
    PenningTrap trap_Euler = trap;

    // Read initial conditions from the trap (first particle) -----
    const auto& P  = trap.particles[0];
    const double q = P.q;
    const double m = P.m;

    const double x0  = P.r(0);
    const double y0  = P.r(1);        
    const double z0  = P.r(2);

    const double vx0 = P.v(0);        
    const double vy0 = P.v(1);        
    const double vz0 = P.v(2);        

    const double B0  = trap.B0;
    const double V0  = trap.V0;
    const double d   = trap.d;

    // Derived frequencies -----
    const double omega0 = q * B0 / m;
    const double omegaz2 = 2.0 * q * V0 / (m * d * d);
    const double omegaz  = std::sqrt(std::max(0.0, omegaz2));
    const double disc = std::max(0.0, 1.0 - 2.0 * omegaz2 / (omega0 * omega0));
    const double omega_plus  = 0.5 * omega0 * (1.0 + std::sqrt(disc));
    const double omega_minus = 0.5 * omega0 * (1.0 - std::sqrt(disc));

    // A± 
    const double denom = (omega_minus - omega_plus);
    const double Aplus  = (vy0 + omega_minus * x0) / denom;
    const double Aminus = -(vy0 + omega_plus  * x0) / denom;

    // Time step 
    const double h = T_us / static_cast<double>(Nstep);

    // Output file 
    std::string filename = "../scripts/r_error_" + std::to_string(Nstep) + ".txt";
    std::ofstream out(filename);
    out << std::scientific << std::setprecision(12);

    // Precompute imaginary unit 
    const std::complex<double> I(0.0, 1.0);

    // Loop over steps
    for (int k = 0; k <= Nstep; ++k) {
        const double t = k * h;

        // z(t) = z0 cos(ωz t)
        const double z = z0 * std::cos(omegaz * t);

        // f(t) = x + i y = A+ e^{-i ω+ t} + A- e^{-i ω- t}
        std::complex<double> f =
            Aplus  * std::exp(-I * omega_plus  * t) +
            Aminus * std::exp(-I * omega_minus * t);

        const double x = std::real(f);
        const double y = std::imag(f);

        // r_exact vector magnitude
        const double r_abs = std::sqrt(x*x + y*y + z*z);

        // RK4 value
        trap_RK4.evolve_RK4(h, 0);
        arma::vec r_diff_RK4(3, arma::fill::zeros);
        r_diff_RK4(0) = x - trap_RK4.particles[0].r(0);
        r_diff_RK4(1) = y - trap_RK4.particles[0].r(1);
        r_diff_RK4(2) = z - trap_RK4.particles[0].r(2);
        const double error_RK4 = arma::norm(r_diff_RK4) / r_abs;  // |r_exact - r_RK4|/|r_exact|
        
        // Euler value
        trap_Euler.evolve_forward_Euler(h, 0);
        arma::vec r_diff_Euler(3, arma::fill::zeros);
        r_diff_Euler(0) = x - trap_Euler.particles[0].r(0);
        r_diff_Euler(1) = y - trap_Euler.particles[0].r(1);
        r_diff_Euler(2) = z - trap_Euler.particles[0].r(2);
        const double error_Euler = arma::norm(r_diff_Euler) / r_abs;  // |r_exact - r_Euler|/|r_exact|
        
        // write: time  |r_exact| rx ry rz for each method (10 columns)
        out << t << " " << error_RK4 << " " << error_Euler << "\n";
    }

    out.close();
}