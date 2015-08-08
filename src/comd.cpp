#include "comd.h"
#include <iostream>

CoMD::CoMD(int _ncheck, Particles *part, Units *units){
  ncheck = _ncheck;
  f = (double *) malloc(sizeof(double) * part->N);
  for (int i = 0; i < part->N; i++)
    f[i] = 0.0;
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
    f[i] = 0.0;

  /* Interaction of particles in the same cell */
  for (int k = 0; k < cells->ncells; k++) {
    Cell *cell = cells->list + k;
    for (int i = 0; i < cell->natoms - 1; i++) {
      int ii = cell->idxlist[i];
      for (int l=0; l < 3; l++) {
	p[l] = part->v[3*ii+l];
	x[l] = part->x[3*ii+l];
      }
      for (int j = i+1; j < cell->natoms; j++) {
	int jj = cell->idxlist[j];
	if (part->spin[ii] != part->spin[jj]) continue;
	if (part->isospin[ii] != part->isospin[jj]) continue;
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
	    occ *= gamma(dp[l]/(sqrt(2)* part->sigma_p));
	    occ *= gamma(dx[l]/(sqrt(2)* part->sigma_r));
	  
	  }
	  f[ii] += occ;	
	  f[jj] += occ;
	}
      }
    }
  }
  
  /* Interaction of particles in different cells */
  for (int k = 0; k < cells->npairs; k++) {
    Cell *c1 = cells->list + cells->pairs[2*k+0];
    Cell *c2 = cells->list + cells->pairs[2*k+1];
    
    for (int i = 0; i < c1->natoms; i++) {
      int ii = c1->idxlist[i];
      for (int l=0; l < 3; l++) {
	p[l] = part->v[3*ii+l];
	x[l] = part->x[3*ii+l];
      }
      
      for (int j = 0; j < c2->natoms; j++) {
	int jj = c2->idxlist[j];
	if (part->spin[ii] != part->spin[jj]) continue;
	if (part->isospin[ii] != part->isospin[jj]) continue;

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
	  int occ = 1.0;
	  for (int l = 0; l < 3; l++) {
	    occ *= gamma(dp[l]/(sqrt(2)* part->sigma_p));
	    occ *= gamma(dx[l]/(sqrt(2)* part->sigma_r));
	  }
	  f[ii] += occ;	
	  f[jj] += occ;
	}
      }
    }
  }
}

