#include "potential.h"

Potential::Potential(Particles *part) {
  /* The default constructor has 2 types of particles and sets the
   values as in CoMD. We pass the particles to also know the size in
   q-space. */

  int n = part->ntypes;

  t0 = -356.0;
  t3 = 303.0;
  r0 = 0.16;
  asym = 32.0;
  cs = -0.33;
  e = 1.2;
  sigma_r = part->sigma_r;
  u = 7.0/6.0;

  rcut = (double **) malloc((n+1)*sizeof(double *));
  for (int i = 0; i < n+1; i++) {
    rcut[i] = (double *) malloc((n+1)*sizeof(double));
  }

  phicut = (double **) malloc((n+1)*sizeof(double *));
  for (int i = 0; i < n+1; i++) {
    phicut[i] = (double *) malloc((n+1)*sizeof(double));
    
  }

  rcut[1][1] = 5.0;
  rcut[2][1] = 5.0;
  rcut[1][2] = 5.0;
  rcut[2][2] = 5.0;

  phicut[1][1] = 0.0;
  phicut[2][1] = 0.0;
  phicut[1][2] = 0.0;
  phicut[2][2] = 0.0;

  
  dphi(rcut[1][1], 1, 1, &(phicut[1][1]));
  dphi(rcut[2][1], 2, 1, &(phicut[2][1]));
  dphi(rcut[1][2], 1, 2, &(phicut[1][2]));
  dphi(rcut[2][2], 2, 2, &(phicut[2][2]));
}

double Potential::rho(double r) {
  return 1.0/pow(2*sqrt(M_PI)*sigma_r,3)* exp(-r*r/(4*sigma_r*sigma_r));
}

double Potential::dphi(double r, int t1, int t2, double *pe) {
  /* Need to work out the expression for dphi */

  double vol = t0/(4*r0) * rho(r)/(sigma_r*sigma_r);
  double three = t3 * u/(u +1) * pow(rho(r)/r0, u)/(sigma_r*sigma_r);
  double sym = asym/4 * (rho(r)/r0) * 1.0/(sigma_r*sigma_r);
  double surf = cs/(4*pow(sigma_r,4)) * (rho(r)/r0) * (r*r/(2*sigma_r*sigma_r) - 3);
  double coul = -e*e/(r*r) * ((M_PI * sigma_r * sigma_r * 8) * rho(r) - erf(r/(2*sigma_r))/r);
  
  
  double evol = t0/2 * rho(r)/r0;
  double ethree = 2*t3/(u +1) * pow(rho(r)/r0, u);
  double esym = asym/2 * rho(r)/r0;
  double esurf = cs/(4*pow(sigma_r,4)) * rho(r)/r0 * (r*r - 2*sigma_r * sigma_r);
  double ecoul = e*e/r * erf(r/(2*sigma_r));



  if (t1 == 1 and t2 == 1) {
    *pe = evol + ethree + esym + esurf;
    return vol + three + sym + surf;
  }
  else if (t1 == 2 and t2 == 2) {
    *pe = evol + ethree + esym + esurf + ecoul;
    return vol + three + sym + surf + coul;
  } else {
    *pe = evol + ethree - esym + esurf;
    return vol + three - sym + surf;
  }


  // double invdr = 1.0/dr;
  // double invdr2 = invdr * invdr;
  // double invdr6 = invdr2 * invdr2 * invdr2;
  // double invdr12 = invdr6 * invdr6;
  // if (t1 != t2) {
  //   invdr12 *= 0;
  //   invdr6 *= 0;
  // }
  // *pe = 4.0 * (invdr12 - invdr6) - phicut[t1][t2];
  // return 24.0 * invdr2 * (2.0 * invdr12 - invdr6);
}

void Potential::write_table(int npoints, double rmin, double rmax) {
  std::ofstream file;
  double r, etot;
  file.open("table.dat");
  file << "#r, Vvol, Vthree, Vsym, Vsurf, Fvol, Fthree, Fsym, Fsurf" << std::endl;
  for (int i = 0; i < npoints; i++){
    r = (rmax - rmin) * (float)(i + 1)/npoints + rmin;
    double vol = t0/(4*r0) * rho(r)/(sigma_r*sigma_r);
    double three = t3 * u/(u +1) * pow(rho(r)/r0, u)/(sigma_r*sigma_r);
    double sym = asym/4 * (rho(r)/r0) * 1.0/(sigma_r*sigma_r);
    double surf = cs/(4*pow(sigma_r,4)) * (rho(r)/r0) * (r*r/(2*sigma_r*sigma_r) - 3);

    
    double evol = t0/2 * rho(r)/r0;
    double ethree = 2*t3/(u +1) * pow(rho(r)/r0, u);
    double esym = asym/2 * rho(r)/r0;
    double esurf = cs/(4*pow(sigma_r,4)) * rho(r)/r0 * (r*r - 2*sigma_r * sigma_r);


    file << r << " ";
    file << evol << " " << ethree << " " << esym << " " << esurf << " ";
    file << r * vol << " " << r * three << " " << r * sym << " " << r * surf << std::endl;
  }
  return;
}
    
    
    
    
