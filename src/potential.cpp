#include "potential.h"

Potential::Potential(Particles *part) {
  /* The default constructor has 2 types of particles and sets the
   values as in CoMD. We pass the particles to also know the size in
   q-space. */

  r0 = ;
  p0 = ;
  v0 = ;
  hbarra = ;
  D = ;
  
  
  int n = part->ntypes;



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


double Potential::dphi(double r, double p, int t1, int t2, double *pe) {

  return -(r/(r0^2))*v0*((hbarra/(r0*p0))^D)*pow(e,(-p*p/(2*p0^2)))*pow(e,(-r*r/(2*r0^2)));
}

double Potential::dpsi(double r, double p, int t1, int t2, double *pe) {

  return -(p/(p0^2))*v0*((hbarra/(r0*p0))^D)*pow(e,(-p*p/(2*p0^2)))*pow(e,(-r*r/(2*r0^2)));
}




    
    
    
    
