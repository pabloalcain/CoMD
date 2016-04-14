#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "particles.h"

class Integrator {
    /* This is a verlet integrator */
protected:
  friend class System;
  double dt;
  int nstep;
public:
  Integrator(double _dt);
  void first_step(Particles *part); 
  virtual void final_step(Particles *part);
  void kinetic(Particles *part);
};

class Cooldown: public Integrator {
  double alpha;
  int nfreq;

public:
  Cooldown(double _dt, double _alpha, int _nfreq);
  virtual void final_step(Particles *part);
};
#endif
