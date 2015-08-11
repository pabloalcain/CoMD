#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "particles.h"

class Integrator {
    /* This is a verlet integrator */
  friend class System;
  double dt;
  int nfreq;
 public:
  Integrator(double _dt, int nfreq);
  void first_step(Particles *part); 
  void final_step(Particles *part);
  void scale(Particles *part, double alpha);
};       
#endif
