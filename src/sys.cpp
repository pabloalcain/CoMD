#include "sys.h"

System::System(Box *_box, Particles *_part, Potential *_pot,
	       Integrator *_integ, Dump *_dump, Thermo *_thermo){
  box = _box;
  part = _part;
  pot = _pot;
  integ = _integ;
  dump = _dump;
  thermo = _thermo;
  units = new Units();
  part->sigma_p *= units->hbar;
  comd = new CoMD(1, part, units);
  neighbor = new Neighbor(part, pot, 10.0, 1000);
  cells = new CellList(neighbor->maxcutoff, part, neighbor->maxcutoff, box);
}


void System::forces_neigh() {
  double x[3];
  double dx[3];
  double dr, dphi, pe;
  int t1, t2;

  part->pe = 0;
  for (int ii = 0; ii < 3*part->N; ii++)
    part->f[ii] = 0;
  
  for (int ii = 0; ii < part->N; ii++) {
    x[0] = part->x[3*ii+0];
    x[1] = part->x[3*ii+1];
    x[2] = part->x[3*ii+2];
    t1 = part->isospin[ii]?2:1;

    for (int j = 0; j< neighbor->num[ii]; j++) {
      int jj  = neighbor->list[ii*part->N + j];
      t2 = part->isospin[jj]?2:1;
      for (int l = 0; l < 3; l++) {
	dx[l] = box->pbc(x[l] - part->x[3*jj + l], l);
      }
      dr = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (dr < pot->rcut[t1][t2] * pot->rcut[t1][t2]) {
	dr = sqrt(dr);
	dphi = pot->dphi(dr, t1, t2, &pe);
	part->pe += pe;
	for (int l = 0; l < 3; l++) {
	  /* Newton's 3rd law */
	  part->f[3*ii + l] += dphi * dx[l];
	  part->f[3*jj + l] -= dphi * dx[l];
	}
      }
    }
  }
}

void System::forces_all() {
  double *x;
  double dx[3];
  double dr, dphi;
  double pe;
  int t1, t2;
  
  part->pe = 0;
  for (int ii = 0; ii < 3*part->N; ii++)
      part->f[ii] = 0;
  
  for (int ii = 0; ii < part->N-1; ii++) {
    x = part->x + 3*ii;
    t1 = part->isospin[ii]?2:1;
    
    for (int jj = ii+1; jj < part->N; jj++) {
      t2 = part->isospin[jj]?2:1;
      for (int l = 0; l < 3; l++) {
	dx[l] = box->pbc(x[l] - part->x[3*jj + l], l);
      }
      dr = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      

      if (dr < pot->rcut[t1][t2] * pot->rcut[t1][t2]) {
	  dr = sqrt(dr);
	  dphi = pot->dphi(dr, t1, t2, &pe);
	  part->pe += pe;
	  for (int l = 0; l < 3; l++) {
	      /* Newton's 3rd law */
	      part->f[3*ii + l] += dphi * dx[l];
	      part->f[3*jj + l] -= dphi * dx[l];
	  }
      }
    }
  }
}

void System::run(int nsteps) {
  cells->update(part, box);
  neighbor->update(cells, part, pot, box);
  integ->kinetic(part);
  forces_neigh();
    
  for (int i = 0; i < nsteps; i++) {
    if (i % comd->ncheck == 0) comd->check_occupation(part, cells, box);
    if (i % dump->nfreq == 0) dump->write(i, part, box);
    if (i % thermo->nfreq == 0) thermo->write(i, part);

    integ->first_step(part);
    if (i % neighbor->nfreq == 0) neighbor->update(cells, part, pot, box);
    forces_neigh();
    integ->final_step(part);
  }
}

