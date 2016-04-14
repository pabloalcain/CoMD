#include "integrator.h"
#include <iostream>

Integrator::Integrator(double _dt) {
    dt = _dt;
    nstep = 0;
}

void Integrator::first_step(Particles *part) {
  double dtmf =  0.5 * dt / part->mass;

  nstep++;
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

void Integrator::kinetic(Particles *part) {
    part->ke = 0;
    for (int i = 0; i < 3*part->N; i++) {
	part->ke += (part->v[i] * part->v[i]);
    }
    part->ke *= 0.5 * part->mass;
}  

Cooldown::Cooldown(double _dt, double _alpha, int _nfreq)
  : Integrator(_dt) {
    nfreq = _nfreq;
    alpha = _alpha;
}

void Cooldown::final_step(Particles *part) {
  Integrator::final_step(part);
  if (nstep % nfreq == 0) {
    for (int i = 0; i < 3*part->N; i++) {
      part->v[i] *= alpha;
    }
  }
}
