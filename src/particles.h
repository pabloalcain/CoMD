#ifndef PARTICLES_H
#define PARTICLES_H

#include "box.h"

#include <math.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

class Particles    
{
  /* This is the particles structure. We create a SoA instead of a
     AoS. The variables in here are those that are needed to
     characterize a wave packet in CoMD: 
     http://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.024612

     The isospin is false for neutrons, true for protons.
  */
  friend class Integrator;
  friend class Cooldown;
  friend class System;
  friend class CellList;
  friend class Dump;
  friend class CoMD;
  friend class Neighbor;
  friend class Potential;
  
  int N;
  int ntypes;
  double *x, *v, *f, *occ;
  bool *isospin, *spin;
  double sigma_r, sigma_p;
  double mass;
  void allocate();
  

public:
  Particles(int nprot, int nneut, Box *box);
  Particles(const std::string& fname, Box *box);
  double pe, ke;
};
#endif

