#include "sys.h"

System::System(Box *_box, Particles *_part, Potential *_pot,
	       Integrator *_integ, Dump *_dump, Thermo *_thermo){
  box = _box;
  part = _part;
  pot = _pot;
  integ = _integ;
  dump = _dump;
  thermo = _thermo;
  cells = new CellList(pot->rcut, part, pot, box);
}

void System::forces() {
  double x[3];
  double dx[3];
  double dr, dphi, pe;

  part->pe = 0;

  for (int ii = 0; ii < 3*part->N; ii++)
      part->f[ii] = 0;


  /* Interaction of particles in the same cell */
  for (int k = 0; k < cells->ncells; k++) {
    Cell *cell = cells->list + k;
    for (int i = 0; i < cell->natoms - 1; i++) {
      int ii = cell->idxlist[i];
      for (int l = 0; l < 3; l++) {
	x[l] = part->x[3*ii + l];
      }
      
      for (int j = i+1; j < cell->natoms; j++) {
	int jj = cell->idxlist[j];
	for (int l = 0; l < 3; l++) {
	  dx[l] = x[l] - part->x[3*jj + l];
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

  /* Interaction of particles in different cells */
  for (int k = 0; k < cells->npairs; k++) {
    Cell *c1 = cells->list + cells->pairs[2*k+0];
    Cell *c2 = cells->list + cells->pairs[2*k+1];
    for (int i = 0; i < c1->natoms; i++) {
      int ii = c1->idxlist[i];
      for (int l = 0; l < 3; l++) {
	x[l] = part->x[3*ii + l];
      }
      
      for (int j = 0; j < c2->natoms; j++) {
	int jj = c2->idxlist[j];
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
}  

void System::forces_all() {
  double x[3];
  double dx[3];
  double dr, dphi;
  double pe;

  part->pe = 0;
  for (int ii = 0; ii < 3*part->N; ii++)
      part->f[ii] = 0;
  
  for (int ii = 0; ii < part->N-1; ii++) {
    for (int l = 0; l < 3; l++) {
      x[l] = part->x[3*ii + l];
    }
    
    for (int jj = ii+1; jj < part->N; jj++) {
      for (int l = 0; l < 3; l++) {
	dx[l] = x[l] - part->x[3*jj + l];
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
    
  for (int i = 0; i < nsteps; i++) {
    if (i % dump->nfreq == 0) {
	dump->write(i, part, box);
    }
    if (i % thermo->nfreq == 0) {
	thermo->write(i, part);
    }
    integ->first_step(part);
    forces_all();
    integ->final_step(part);
    cells->update(part, box);
    
  }
}

