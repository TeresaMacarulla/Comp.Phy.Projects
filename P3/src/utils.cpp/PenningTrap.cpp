#include "PenningTrap.hpp"

// Definitions of constructors
PenningTrap::PenningTrap(double B0, double V0, double d, std::vector<Particle>)
{
  // ...
}

void PenningTrap::add_particle(Particle p_in)
{
}

arma::vec PenningTrap::external_E_field(arma::vec r)
{
  const double coeff = V0 / (d * d);
  arma::vec E(3, arma::fill::none);
  E(0) = coeff * r(0);
  E(1) = coeff * r(1);
  E(2) = -2.0 * coeff * r(2);
  return E;
}

arma::vec PenningTrap::external_B_field(arma::vec r);
{
  arma::vec B(3, arma::fill::none);
  B(0) = 0.0;
  B(1) = 0.0;
  B(2) = B0;
  return B;
}

arma::vec PenningTrap::force_particle(int i, int j);
{

}

arma::vec PenningTrap::total_force_external(arma::vec r, arma::vec v);
{
  // Obtain E and B fields from previous functions 
  arma::vec E = external_E_field(r);
  arma::vec B = external_B_field(r);

  // Calculate Lorentz force: F = q (E + v × B)
  arma::vec F = Particla[1].q * (E + arma::cross(v, B));

  return F;
}

arma::vec PenningTrap::total_force_particles(int i);
{
}

arma::vec PenningTrap::total_force(int i);
{
  
}

void PenningTrap::evolve_RK4(double h);
{
  // Number of particles in the Penning Trap: n
  const std::size_t n = Particle.size();

  // Previous condition 
  std::vector<arma::vec> r0(n), v0(n);
  for (std::size_t i = 0; i < n; ++i) {
      r0[i] = Particle[i].r;
      v0[i] = Particle[i].v;
  }

  // k for r and v (each one k(i) is arma::vec with 3 components)
  std::vector<arma::vec> kr1(n), kr2(n), kr3(n), kr4(n);
  std::vector<arma::vec> kv1(n), kv2(n), kv3(n), kv4(n);

  // ---- k1 ----
  for (std::size_t i = 0; i < n; ++i) {
    kr1[i] = h * v0[i];
    kv1[i] = h * (total_force(r0[i],v0[i]) / Particle[i].m);
  }

  // ---- k2 ----
  std::vector<arma::vec> r12(n), v12(n);
  for (std::size_t i = 0; i < n; ++i) {
      r12[i] = r0[i] + 0.5 * kr1[i];
      v12[i] = v0[i] + 0.5 * kv1[i];
      kr2[i] = h * v12[i];
      kv2[i] = h * (total_force(r12[i],v12[i]) / Particle[i].m);
  }

  // ---- k3 ----
  std::vector<arma::vec> r32(n), v32(n);
  for (std::size_t i = 0; i < n; ++i) {
    r32[i] = r0[i] + 0.5 * kr2[i];
    v32[i] = v0[i] + 0.5 * kv2[i];
    kr3[i] = h * v32[i];
    kv3[i] = h * (total_force(r32[i],v32[i]) / Particle[i].m);
  }

  // ---- k4 ----
  std::vector<arma::vec> r4(n), v4(n);
  for (std::size_t i = 0; i < n; ++i) {
    r4[i] = r0[i] + kr3[i];
    v4[i] = v0[i] + kv3[i];
    kr4[i] = h * v4[i];
    kv4[i] = h * (total_force(r4[i],v4[i]) / Particle[i].m);
  }

  // ---- final actualization ----
  for (std::size_t i = 0; i < n; ++i) {
    Particle[i].r = r0[i] + (kr1[i] + 2.0*kr2[i] + 2.0*kr3[i] + kr4[i]) / 6.0;
    Particle[i].v = v0[i] + (kv1[i] + 2.0*kv2[i] + 2.0*kv3[i] + kv4[i]) / 6.0;
  }
}

void PenningTrap::evolve_forward_Euler(double h);
{
  const std::size_t n = Particle.size();

  // Previous velocity condition 
  std::vector<arma::vec> v_old(n);
  for (std::size_t i = 0; i < n; ++i){
    v_old[i] = Particle[i].v;
    r_old[i] = Particle[i].r;
  } 

  for (std::size_t i = 0; i < n; ++i) {
    Particle[i].v = v_old[i] + h * (total_force(r_old[i],v_old[i]) / Particle[i].m);
    Particle[i].r = r_old[i] + h * v_old[i];
  }
}