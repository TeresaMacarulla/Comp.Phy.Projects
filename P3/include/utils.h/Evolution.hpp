#ifndef __EVOLUTION_HPP__
#define __EVOLUTION_HPP__

#include "PenningTrap.hpp"
#include "PenningTrap_Vt.hpp"
#include <cstddef> // para std::size_t
#include <iomanip>
#include <iostream> 
#include <armadillo> 
#include <fstream>
#include <complex>

// Function declarations

// Simulate a single particle and track z-direction movement
void evolution_rz(PenningTrap& trap, double T_us, int Nstep);

// Simulate two particles and track rx and ry movement 
void evolution_two(PenningTrap& trap, double T_us, int Nstep, int inter = 1);

// Applies analytical solution, RK4 and Euler method. Writes the relative error of Euler and RK4 to a file 
void evolution_RK4_Euler(PenningTrap& trap, double T_us, int Nstep);

// Compare analytical solution, RK4 and Euler for single particle and runs different nsteps (4000, 8000, 16000, 32000)
void evolution_nsteps(int n, PenningTrap& trap);

// Fill trap with particles with randomly generated initial positions and velocities
void fill_trap(PenningTrap_Vt& trap_Vt, double omegav_in, double f_in, int N);

// Simulates many particle system 
void evolution_many(PenningTrap_Vt& trap, double T_us, int Nstep, int inter);

#endif
