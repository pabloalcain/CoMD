#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "math.h"
#include <iostream>

#include "particles.h"


class Potential
{
  /* Basic CoMD potential structure, with all the parameters */
  double t0, t3;
  double r0;
  double asym;
  double cs;
  double e;
  
  double sigma_r;
  double u;
  
  double rho(double r);
    
public:
  double **rcut, **phicut;
  Potential(Particles *part);
  double dphi(double r, int t1, int t2, double *pe);
  void write_table(int npoints, double rmin, double rmax);

};

#endif
