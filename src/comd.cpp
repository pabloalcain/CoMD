#include "comd.h"

CoMD::CoMD(int _ncheck, Particles *part, Units *units){
  ncheck = _ncheck;
  create_lut(1000, part->sigma_r, part->sigma_p, units->hbar);
}

void CoMD::create_lut(int npoints, double sr, double sp, double hbar) {
  double u, du;

  lut_gamma = (double *) malloc(sizeof(double) * (npoints));
  lut_rmax = 3/sqrt(2);
  lut_invrmax = 1.0/lut_rmax;
  lut_npoints = npoints;
  for (int i = 0; i < npoints; i++) {
    u = (double)i * lut_rmax/npoints;
    du = sqrt(2*M_PI*hbar/(2 * sr * sp));
    lut_gamma[i] = (erf(u + du) - erf(u - du))/2;
  }
}

double CoMD::gamma(double u) {
  /* This gamma function comes from the integration of the f given in
     the paper */
  int idx = (int) (lut_npoints * u * lut_invrmax);
  return lut_gamma[idx];
}

void CoMD::check_occupation(Particles *part, CellList *cells, Box *box) {

  double p[3], dp[3], sqp;
  double x[3], dx[3], sqx;
  
  for (int i = 0; i < part->N; i++)
    part->occ[i] = lut_gamma[0];

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
}

