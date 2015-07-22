#include "integrator.h"

Integrator::Integrator(double _dt) {
    dt = _dt;
}

void Integrator::first_step(Particles *part) {
    double dtmf =  0.5 * dt / part->mass;
    
    for (int i = 0; i < 3*part->N; i++) {
	part->v[i] += dtmf * part->f[i];
	part->x[i] += dt * part->v[i];
    }
}

void Integrator::final_step(Particles *part) {
    double dtmf =  0.5 * dt / part->mass;
    
    for (int i = 0; i < 3*part->N; i++) {
	part->v[i] += dtmf * part->f[i];
    }
}
