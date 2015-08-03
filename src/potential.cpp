#include "potential.h"

Potential::Potential() {
  /* Dummy values so far */
  t0 = 1.0;
  r0 = 1.0;
  t3 = 1.0;
  asym = 1.0;
  cs = 1.0;
  e = 1.0;
  rcut = 1.0;
}

double Potential::dphi(double dr) {
  /* Need to work out the expression for dphi */
  double invdr = 1.0/dr;
  double invdr2 = invdr * invdr;
  double invdr4 = invdr2 * invdr2;
  double invdr6 = invdr4 * invdr2;
  double invdr12 = invdr6 * invdr6;
  return 24 * invdr2 * (2.0*invdr12 - invdr6);
}
	
