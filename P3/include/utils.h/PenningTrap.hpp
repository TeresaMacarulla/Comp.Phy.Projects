#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <armadillo>
#include <vector>
#include "Particle.hpp"

class PenningTrap 
{

public:

  PenningTrap(double B0, double V0, double d, std::vector<Particle>);

  double electric_field(arguments);
  double magnetic_field(arguments);
  double force(arguments);

}; 

#endif