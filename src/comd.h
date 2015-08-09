#ifndef COMD_H
#define COMD_H

#include "cell.h"
#include "particles.h"
#include "box.h"
#include "units.h"

#include <cmath>

class CoMD
{
  /* CoMD class, the key idea of the development. We follow loosely
     prescriptions by Bonasera et al., PRC64 (024612)
     http://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.024612.
     Some changes are: we won't add Pauli blocking, and not necessarily
     check in every step whether f gets larger than 1.
     
     We might change slightly the algorithm, maybe aiming towards a
     minimization of f?


  */

  double gamma(double u);
  double sigma_r, sigma_p;
  void create_lut(int npoints, double sigma_r, double sigma_p, double hbar);
  double *lut_gamma;
  double lut_rmax, lut_invrmax, lut_npoints;
  
 public:
  int ncheck;
  
  CoMD(int ncheck, Particles *part, Units *units);
  void check_occupation(Particles *part, CellList *cells, Box *box);
  void change_momentum(Particles *part, CellList *cells);
};

#endif
