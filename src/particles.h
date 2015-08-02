#ifndef PARTICLES_H
#define PARTICLES_H

#include "box.h"
#include <math.h>
#include <stdlib.h>
class Particles    
{
    /* This is the particles structure. We create a SoA instead of a
       AoS. The variables in here are those that are needed to
       characterize a wave packet in CoMD: 
       http://journals.aps.org/prc/abstract/10.1103/PhysRevC.64.024612

       The isospin is false for neutrons, true for protons.
    */
    friend class Integrator;
    friend class System;
    friend class CellList;
    friend class Dump;
    int N;
    double *x, *v, *f;
    bool *isospin, *spin;
    double sig_r, sig_p;
    double mass;

public:
    Particles(int nprot, int nneut, Box box);
    void pbc(Box box);
};
#endif

