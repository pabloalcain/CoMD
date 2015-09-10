#include "comd.h"

CoMD::CoMD(int _ncheck, Particles *part, Units *units){
  ncheck = _ncheck;

  create_lut(1000, part->sigma_r, part->sigma_p, units->hbar);
}

void CoMD::random_rotation(double u1, double u2, double u3, double *i, double *vector) {
  /* Generates a random rotation matrix following Arvo prescription */
  double cos1 = cos(u1);
  double cos2 = cos(u2);
  double sin1 = sin(u1);
  double sin2 = sin(u2);

  /* rotation matrix */
  /* those we don't set are zeros */
  double r11 = cos1;
  double r12 = sin1;
  double r21 = -sin1;
  double r22 = cos1;
  /* r33 is 1, so we avoid */

  /* householder matrix */
  double h11 = 1 - 2 * u3 * cos2 * cos2;
  double h12 = -2 * u3 * cos2 * sin2;
  double h13 = -2 * cos2 * sqrt(u3 * (1 - u3));
  double h21 = h12;
  double h22 = 1 - 2 * u3 * sin2 * sin2;
  double h23 = -2 * sin2 * sqrt(u3 * (1 - u3));
  double h31 = h13;
  double h32 = h23;
  double h33 = 2 * u3 - 1;

  /* result matrix */
  double m11 = h11 * r11 + h12 * r21;
  double m12 = h11 * r12 + h12 * r22;
  double m13 = h13;
  double m21 = h21 * r11 + h22 * r21;
  double m22 = h21 * r12 + h22 * r22;
  double m23 = h23;
  double m31 = h31 * r11 + h32 * r21;
  double m32 = h31 * r12 + h32 * r22;
  double m33 = h33;

  /* rotate input vector */
  vector[0] = i[0]*m11 + i[1] * m12 + i[2] * m13;
  vector[1] = i[0]*m21 + i[1] * m22 + i[2] * m23;
  vector[2] = i[0]*m31 + i[1] * m32 + i[2] * m33;
}


void CoMD::create_lut(int npoints, double sr, double sp, double hbar) {
  double u, du;

  lut_gamma = (double *) malloc(sizeof(double) * (npoints));
  lut_rmax = 3/sqrt(2);
  lut_invrmax = 1.0/lut_rmax;
  lut_npoints = npoints;
  for (int i = 0; i < npoints; i++) {
    u = (double)i * lut_rmax/npoints;
    du = 4.5*sqrt(M_PI*hbar/(2 * sr * sp));
    lut_gamma[i] = (erf(u + du) - erf(u - du))/2;
  }
}

double CoMD::gamma(double u) {
  /* This gamma function comes from the integration of the f given in
     the paper */
  int idx = (int) (lut_npoints * u * lut_invrmax);
  return lut_gamma[idx];
}

double CoMD::occup_one(int ii, Particles *part, Neighbor *neighbor, Box *box) {
  double dx[3], dp[3], p[3], x[3];
  double sqp, sqx;
  
  for (int l = 0; l < 3; l++) {
    x[l] = part->x[3*ii + l];
    p[l] = part->v[3*ii + l];
  }
  double occupation = lut_gamma[0]*lut_gamma[0]*lut_gamma[0];
  occupation *= occupation;
  for (int j = 0; j< neighbor->num[ii]; j++) {
    int jj  = neighbor->list[ii*part->N + j];
    if (part->isospin[ii] != part->isospin[jj]) continue;
    if (part->spin[ii] != part->spin[jj]) continue;
    sqp = 0.0;
    sqx = 0.0;
    for (int l = 0; l < 3; l++) {
      dp[l] = part->mass * (p[l] - part->v[3*jj + l]);
      dx[l] = box->pbc(x[l] - part->x[3*jj + l], l);
      sqp += dp[l] * dp[l];
      sqx += dx[l] * dx[l];
    }
    /* We check whether the particles are inside the 3sigma
       ellipsoid */
    double u = sqp/(3 * part->sigma_p * 3 * part->sigma_p);
    u += sqx/(3 * part->sigma_r * 3 * part->sigma_r);
    if (u < 1) {
      double occ = 1.0;
      for (int l = 0; l < 3; l++) {
        occ *= gamma(fabs(dp[l])/(sqrt(2)* part->sigma_p));
        occ *= gamma(fabs(dx[l])/(sqrt(2)* part->sigma_r));
      }
      occupation += occ;	
    }
  }
  return occupation;
}  

void CoMD::check_occupation(Particles *part, Neighbor *neighbor, Box *box) {

  double p[3], dp[3], sqp;
  double x[3], dx[3], sqx;
  double occupation = lut_gamma[0]*lut_gamma[0]*lut_gamma[0];
  occupation *= occupation;

  for (int i = 0; i < part->N; i++)
    part->occ[i] = occupation;

  for (int ii = 0; ii < part->N; ii++) {
    x[0] = part->x[3*ii+0];
    x[1] = part->x[3*ii+1];
    x[2] = part->x[3*ii+2];
    
    for (int j = 0; j< neighbor->num[ii]; j++) {
      int jj  = neighbor->list[ii*part->N + j];
      if (part->isospin[ii] != part->isospin[jj]) continue;
      if (part->spin[ii] != part->spin[jj]) continue;
      sqp = 0.0;
      sqx = 0.0;
      for (int l = 0; l < 3; l++) {
        dp[l] = part->mass * (p[l] - part->v[3*jj + l]);
        dx[l] = box->pbc(x[l] - part->x[3*jj + l], l);
        sqp += dp[l] * dp[l];
        sqx += dx[l] * dx[l];
      }
      /* We check whether the particles are inside the 3sigma
         ellipsoid */
      double u = sqp/(3 * part->sigma_p * 3 * part->sigma_p);
      u += sqx/(3 * part->sigma_r * 3 * part->sigma_r);
      if (u < 1) {
        /* std::cout << "Yo ayudo!" << std::endl; */
        double occ = 1.0;
        for (int l = 0; l < 3; l++) {
          occ *= gamma(fabs(dp[l])/(sqrt(2)* part->sigma_p));
          occ *= gamma(fabs(dx[l])/(sqrt(2)* part->sigma_r));
        }
        part->occ[ii] += occ;	
        part->occ[jj] += occ;
      }
    }
  }
}

void CoMD::change_momentum(Particles *part, Neighbor *neighbor, Box *box) {
  int nneigh;

  for (int ii = 0; ii < part->N; ii++) {
    double occ = part->occ[ii];
    while (occ > 1) {
      std::cout << "Estoy cambiando algo!" << std::endl;
      nneigh = neighbor->num[ii];
      /* Draw a random neighbor */
      std::uniform_int_distribution<int> pick(0,nneigh - 1);
      int myneigh = pick(generator);
      int jj = neighbor->list[ii*part->N + myneigh];
      scatter(part, ii, jj);
      occ = occup_one(ii, part, neighbor, box);
    }
    part->occ[ii] = occ;
  }
}

void CoMD::scatter(Particles *part, int ii, int jj){
  /* We consider every particle to have the same mass */

  double dp[3], dprot[3];
  double u1, u2, u3;
  
  for (int l = 0; l < 3; l++)
    dp[l] = part->v[3*ii + l] - part->v[3*jj + l];
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  u1 = distribution(generator);
  u2 = distribution(generator);
  u3 = distribution(generator);
  random_rotation(u1, u2, u3, dp, dprot);
  for (int l = 0; l < 3; l++) {
    part->v[3*ii + l] -= dp[l] - dprot[l];
    part->v[3*jj + l] += dp[l] - dprot[l];
  }
}

