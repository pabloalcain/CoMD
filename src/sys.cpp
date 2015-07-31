#include "sys.h"

System::System(Box *_box, Particles *_part, Potential *_pot, Integrator *_integ){
  box = _box;
  part = _part;
  pot = _pot;
  integ = _integ;
  cells = new CellList(1.0, part, pot, box);
}

void System::forces() {
  double x[3];
  double dx[3];
  double dr, dphi;

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
	  dx[l] = x[l] - part->x[3*ii + l];
	}
	dr = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
	dr = sqrt(dr);
	dphi = pot->dphi(dr);
	for (int l = 0; l < 3; l++) {
	  /* Newton's 3rd law */
	  part->f[3*ii + l] += dphi * dx[l];
	  part->f[3*jj + l] -= dphi * dx[l];
	}
      }
    }
  }

  /* Interaction of particles in different cells */
  for (int k = 0; k < cells->npairs; k++) {
    Cell *c1 = cells->list + 2*k+0;
    Cell *c2 = cells->list + 2*k+1;
    for (int i = 0; i < c1->natoms; i++) {
      int ii = c1->idxlist[i];
      for (int l = 0; l < 3; l++) {
	x[l] = part->x[3*ii + l];
      }
      
      for (int j = 0; j < c2->natoms; j++) {
	int jj = c2->idxlist[j];
	for (int l = 0; l < 3; l++) {
	  dx[l] = x[l] - part->x[3*ii + l];
	}
	dr = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
	dr = sqrt(dr);
	dphi = pot->dphi(dr);
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
  for (int i = 0; i < nsteps; i++) {
    std::cout << "Step " << i << " out of "<< nsteps << std::endl;
    cells->update(part);
    integ->first_step(part);
    forces();
    integ->final_step(part);
  }
}

