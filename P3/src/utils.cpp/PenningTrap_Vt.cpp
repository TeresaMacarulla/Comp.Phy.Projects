#include "PenningTrap_Vt.hpp"

// Definitions of constructors
PenningTrap_Vt::PenningTrap_Vt(double B0_in, double V0_in, double d_in, double omegav_in, double f_in, std::vector<Particle> particles_in)
: B0(B0_in), V0(V0_in), d(d_in), omegav(omegav_in), f(f_in), particles(particles_in) 
{

}

// Add a particle to the trap
void PenningTrap_Vt::add_particle(Particle p_in)
{
  this->particles.push_back(p_in);
}

// External electric field at point r=(x,y,z)
arma::vec PenningTrap_Vt::external_E_field(arma::vec r, double t)
{
  const double coeff = V0 / (d * d);
  arma::vec E(3, arma::fill::none);
  double sum = (1.0 + f * std::cos(omegav * t));

  E(0) = coeff * r(0) * sum;
  E(1) = coeff * r(1) * sum;
  E(2) = -2.0 * coeff * r(2) * sum;

  if (arma::norm(r) > d) {
    arma::vec E(3, arma::fill::none);
  }

  return E;
}

// External magnetic field at point r=(x,y,z)
arma::vec PenningTrap_Vt::external_B_field(arma::vec r)
{
  arma::vec B(3, arma::fill::none);
  B(0) = 0.0;
  B(1) = 0.0;
  B(2) = B0;

  if (arma::norm(r) > d) {
    arma::vec B(3, arma::fill::none);
  }

  return B;
}

// Force on particle_i from particle_j
arma::vec PenningTrap_Vt::force_particle(int i, int j)
{
  arma::vec r_diff(3, arma::fill::zeros);
  arma::vec F(3, arma::fill::zeros);

  double ke = 1.38935333e5;
  r_diff = particles[i].r- particles[j].r;
  //r_diff(1) = particles[i].r(1)- particles[j].r(1);
  //r_diff(2) = particles[i].r(2)- particles[j].r(2);
  double r_norm = arma::norm(r_diff);  // |r_i - r_j|
  F = ke * particles[i].q * particles[j].q * r_diff / std::pow(r_norm, 3);
  return F;
}

// The total force on particle_i from the external fields
arma::vec PenningTrap_Vt::total_force_external(arma::vec r, arma::vec v, double t)
{
  //arma::vec E(3, arma::fill::zeros);
  //arma::vec B(3, arma::fill::zeros);
  //arma::vec F(3, arma::fill::zeros);

  // Obtain E and B fields from previous functions 
  arma::vec E = external_E_field(r, t);
  arma::vec B = external_B_field(r);

  // Calculate Lorentz force: F = q (E + v × B)
  arma::vec F = particles[0].q * (E + arma::cross(v, B));
  return F;
}

// The total force on particle_i from the other particles
arma::vec PenningTrap_Vt::total_force_particles(int i)
{
  double ke = 1.38935333e5;
  arma::vec F(3, arma::fill::zeros);

  for (int j=0; j<particles.size() ; ++j) {
    if (j == i) {
      continue;   // Don't calculate force on particle_i from particle_i
    }
    F = F + force_particle(i,j);
  }

  return F;
}

// The total force on particle_i from both external fields and other particles
arma::vec PenningTrap_Vt::total_force(arma::vec r, arma::vec v, int i, int inter, double t)
{
  arma::vec F(3, arma::fill::zeros);
  const std::size_t n = particles.size();

  if (n==1 || inter==0){
    F = total_force_external(r, v, t);
  }
  else if (n>1){
    F = total_force_external(r, v, t) + total_force_particles(i);
  }
  
  return F;
}

void PenningTrap_Vt::evolve_RK4(double h, int inter, double t)
{
  // Number of particles in the Penning Trap: n
  const std::size_t n = particles.size();

  // Previous condition 
  std::vector<arma::vec> r0(n), v0(n);
  for (std::size_t i = 0; i < n; ++i) {
      r0[i] = particles[i].r;
      v0[i] = particles[i].v;
  }

  // k for r and v (each one k(i) is arma::vec with 3 components)
  std::vector<arma::vec> kr1(n), kr2(n), kr3(n), kr4(n);
  std::vector<arma::vec> kv1(n), kv2(n), kv3(n), kv4(n);

  // ---- k1 ----
  for (std::size_t i = 0; i < n; ++i) {
    kr1[i] = h * v0[i];
    kv1[i] = h * (total_force(r0[i], v0[i], i, inter, t) / particles[i].m);
  }

  // ---- k2 ----
  std::vector<arma::vec> r12(n), v12(n);
  for (std::size_t i = 0; i < n; ++i) {
      r12[i] = r0[i] + 0.5 * kr1[i];
      v12[i] = v0[i] + 0.5 * kv1[i];
      kr2[i] = h * v12[i];
      kv2[i] = h * (total_force(r12[i], v12[i], i, inter, t) / particles[i].m);
  }

  // ---- k3 ----
  std::vector<arma::vec> r32(n), v32(n);
  for (std::size_t i = 0; i < n; ++i) {
    r32[i] = r0[i] + 0.5 * kr2[i];
    v32[i] = v0[i] + 0.5 * kv2[i];
    kr3[i] = h * v32[i];
    kv3[i] = h * (total_force(r32[i], v32[i], i, inter, t) / particles[i].m);
  }

  // ---- k4 ----
  std::vector<arma::vec> r4(n), v4(n);
  for (std::size_t i = 0; i < n; ++i) {
    r4[i] = r0[i] + kr3[i];
    v4[i] = v0[i] + kv3[i];
    kr4[i] = h * v4[i];
    kv4[i] = h * (total_force(r4[i], v4[i], i, inter, t) / particles[i].m);

    // final
    particles[i].r = r0[i] + (kr1[i] + 2.0*kr2[i] + 2.0*kr3[i] + kr4[i]) / 6.0;
    particles[i].v = v0[i] + (kv1[i] + 2.0*kv2[i] + 2.0*kv3[i] + kv4[i]) / 6.0;
  }

  /*
  // ---- final actualization ----
  for (std::size_t i = 0; i < n; ++i) {
    particles[i].r = r0[i] + (kr1[i] + 2.0*kr2[i] + 2.0*kr3[i] + kr4[i]) / 6.0;
    particles[i].v = v0[i] + (kv1[i] + 2.0*kv2[i] + 2.0*kv3[i] + kv4[i]) / 6.0;
  }*/
}

void PenningTrap_Vt::evolve_forward_Euler(double h, int inter, double t)
{
  const std::size_t n = particles.size();

  // Previous velocity condition 
  std::vector<arma::vec> v_old(n), r_old(n);
  for (std::size_t i = 0; i < n; ++i){
    v_old[i] = particles[i].v;
    r_old[i] = particles[i].r;
  } 

  for (std::size_t i = 0; i < n; ++i) {
    particles[i].v = v_old[i] + h * (total_force(r_old[i], v_old[i], i, inter, t) / particles[i].m);
    particles[i].r = r_old[i] + h * v_old[i];
  }
}

std::size_t PenningTrap_Vt::inside_particles()
{
  std::size_t cnt = 0;
  for (const auto& p : particles) {
    if (arma::norm(p.r) < d) {     // |r| = sqrt(x^2+y^2+z^2)
      ++cnt;
    }
  }
  return cnt;
}