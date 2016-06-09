#include "potential.h"
#include <iostream>

Potential::Potential() {
  /* The default constructor has 2 types of particles and sets the
   values as in CoMD. We pass the particles to also know the size in
   q-space. */

  r0 = 1;
  p0 = 1;
  v0 = 1;
  hbarra = 1;
  D = 1;
  
  /*
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
  */
}


double Potential::dphi(double dr, double dp, double *pe) {
  double num = -(dr/ (pow(r0, 2)) ) * v0 * pow( (hbarra/(r0 * p0)), D);
  double exp = pow( e, -dp * dp / (2 * pow(p0, 2)) ) * pow( e, -dr * dr / (2 * pow(r0, 2)) );
  return num * exp;
}

double Potential::dpsi(double dr, double dp, double *pe) {
  double num = -(dp/ (pow(p0, 2)) ) * v0 * pow( (hbarra/(r0 * p0)), D);
  double exp = pow( e, -dp * dp / (2 * pow(p0, 2)) ) * pow( e, -dr * dr / (2 * pow(r0, 2)) );
  return num * exp;
}


    
    
    
    
