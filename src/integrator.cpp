#include "integrator.h"

Integrator::Integrator(double _dt, int _nfreq) {
    dt = _dt;
    nfreq = _nfreq;
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
    part->ke = 0;
    
    for (int i = 0; i < 3*part->N; i++) {
	part->v[i] += dtmf * part->f[i];
	part->ke += (part->v[i] * part->v[i]);
    }
    part->ke *= 0.5 * part->mass;
}

void Integrator::scale(Particles *part, double alpha) {

  for (int i = 0; i < 3*part->N; i++) {
    part->v[i] *= alpha;
  }
}
