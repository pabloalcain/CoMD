#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "particles.h"

class Integrator {
    /* This is a verlet integrator */
    double dt;
public:
    Integrator(double _dt);
    void first_step(Particles *part); 
    void final_step(Particles *part);
};       
#endif
