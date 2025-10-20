#ifndef __PenningTrap_Vt_hpp__
#define __PenningTrap_Vt_hpp__

#include <armadillo>
#include <vector>
#include <cmath>
#include "Particle.hpp"

class PenningTrap_Vt
{

public:

  double B0;
  double V0;
  double d;
  double omegav;
  double f;
  std::vector<Particle> particles;

  PenningTrap_Vt(double B0_in, double V0_in, double d_in, double omegav_in, double f_in, std::vector<Particle> particles_in);

  // Add a particle to the trap
  void add_particle(Particle p_in);

  // External electric field at point r=(x,y,z)
  arma::vec external_E_field(arma::vec r, double t); 

  // External magnetic field at point r=(x,y,z)
  arma::vec external_B_field(arma::vec r);

  // Force on particle_i from particle_j
  arma::vec force_particle(int i, int j);

  // The total force on particle_i from the external fields
  arma::vec total_force_external(arma::vec r, arma::vec v, double t);

  // The total force on particle_i from the other particles
  arma::vec total_force_particles(int i);

  // The total force on particle_i from both external fields and other particles
  arma::vec total_force(arma::vec r, arma::vec v, int i, int inter, double t);

  // Evolve the system one time step (dt) using Runge-Kutta 4th order
  void evolve_RK4(double h, int inter, double t);

  // Evolve the system one time step (dt) using Forward Euler
  void evolve_forward_Euler(double h, int inter, double t);

  // Count the number of particles inside the trap
  std::size_t inside_particles();

}; 

#endif