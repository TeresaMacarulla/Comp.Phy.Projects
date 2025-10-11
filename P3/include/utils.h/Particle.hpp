#ifndef __Particle_hpp__
#define __Particle_hpp__

#include <armadillo>

class Particle 
{

public:

  Particle(double q, double m, arma::vec v, arma::vec r);

  void some_function(arguments);

}; 

#endif