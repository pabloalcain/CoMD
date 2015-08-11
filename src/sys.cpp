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
  comd = new CoMD(1, part, units);
  neighbor = new Neighbor(part, pot, 5.0, 1);
  cells = new CellList(neighbor->cutoff, part, neighbor->cutoff, box);
}


void System::forces_neigh() {
  double x[3];
  double dx[3];
  double dr, dphi, pe;

  part->pe = 0;
  for (int ii = 0; ii < 3*part->N; ii++)
    part->f[ii] = 0;
  
  for (int ii = 0; ii < part->N; ii++) {
    x[0] = part->x[3*ii+0];
    x[1] = part->x[3*ii+1];
    x[2] = part->x[3*ii+2];
    for (int j = 0; j< neighbor->num[ii]; j++) {
      int jj  = neighbor->list[ii*part->N + j];
      for (int l = 0; l < 3; l++) {
	dx[l] = box->pbc(x[l] - part->x[3*jj + l], l);
      }
      dr = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (dr < pot->rcut * pot->rcut) {
	dr = sqrt(dr);
	dphi = pot->dphi(dr, &pe);
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

  part->pe = 0;
  for (int ii = 0; ii < 3*part->N; ii++)
      part->f[ii] = 0;
  
  for (int ii = 0; ii < part->N-1; ii++) {
    x = part->x + 3*ii;
    
    for (int jj = ii+1; jj < part->N; jj++) {
      for (int l = 0; l < 3; l++) {
	dx[l] = box->pbc(x[l] - part->x[3*jj + l], l);
      }
      dr = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      

      if (dr < pot->rcut * pot->rcut) {
	  dr = sqrt(dr);
	  dphi = pot->dphi(dr, &pe);
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

