#include "potential.h"

Potential::Potential() {
  /* Dummy values so far */
  t0 = 1.0;
  r0 = 1.0;
  t3 = 1.0;
  asym = 1.0;
  cs = 1.0;
  e = 1.0;
  rcut = 3.0;
}

double Potential::dphi(double dr, double *pe) {
  /* Need to work out the expression for dphi */
  double invdr = 1.0/dr;
  double invdr2 = invdr * invdr;
  double invdr6 = invdr2 * invdr2 * invdr2;
  double invdr12 = invdr6 * invdr6;
  *pe = 0*4.0 * (invdr12 - invdr6);
  return 0.0*24.0 * invdr2 *(2.0 * invdr12 - invdr6);
}
