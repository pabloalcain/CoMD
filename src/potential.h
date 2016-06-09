#ifndef POTENTIAL_H
#define POTENTIAL_H
#include "math.h"

class Potential
{
    /* Basic CoMD potential structure, with all the parameters */
    double t0, r0;
    double t3;
    double asym;
    double cs;
    double e;
    double p0;
    double v0;
    double hbarra;
    double D;
    
public:
    double rcut;
    double phicut;
    Potential();
    double dphi(double dr, double dp, double *pe);
    double dpsi(double dr, double dp, double *pe);
};

#endif
