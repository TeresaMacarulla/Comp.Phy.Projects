#include "PenningTrap.hpp"
#include <armadillo>

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
}