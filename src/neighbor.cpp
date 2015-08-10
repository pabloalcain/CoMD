#include "neighbor.h"

Neighbor::Neighbor(Particles *part) {
  list = (int *) malloc(part->N * part->N * sizeof(int));
  num = (int *) malloc(part->N * sizeof(int));
}

void Neighbor::update(CellList *cells, Particles *part, Potential *pot, Box *box) {
  /* Add the atom idx in position "natoms" */
  double x[3];
  double dr;
  int n;

  for (int i = 0; i < part->N; i++) num[i] = 0;
  for (int k = 0; k < cells->ncells; k++) {
    Cell *cell = cells->list + k;
    for (int i = 0; i < cell->natoms - 1; i++) {
      int ii = cell->idxlist[i];
      n = 0;
      x[0] = part->x[3*ii+0];
      x[1] = part->x[3*ii+1];
      x[2] = part->x[3*ii+2];
      
      for (int j = i+1; j < cell->natoms; j++) {
	int jj = cell->idxlist[j];
	dr = 0;
	for (int l = 0; l < 3; l++)
	  dr += box->pbc(x[l] - part->x[3*jj + l], l) * box->pbc(x[l] - part->x[3*jj + l], l);
	if (dr < pot->rcut * pot->rcut){
	  list[ii*part->N + n] = jj;
	  n++;
	}
      }
      num[ii] = n;
    }
  }
    
  /* Interaction of particles in different cells */
  for (int k = 0; k < cells->npairs; k++) {
    Cell *c1 = cells->list + cells->pairs[2*k+0];
    Cell *c2 = cells->list + cells->pairs[2*k+1];
    for (int i = 0; i < c1->natoms; i++) {
      int ii = c1->idxlist[i];
      n = num[ii];
      x[0] = part->x[3*ii+0];
      x[1] = part->x[3*ii+1];
      x[2] = part->x[3*ii+2];
      for (int j = 0; j < c2->natoms; j++) {
	int jj = c2->idxlist[j];
	dr = 0;
	for (int l = 0; l < 3; l++)
	  dr += box->pbc(x[l] - part->x[3*jj + l], l) * box->pbc(x[l] - part->x[3*jj + l], l);

	if (dr < pot->rcut * pot->rcut){
	  list[ii*part->N + n] = jj;
	  n++;
	}
      }
      num[ii] = n;
    }
  }
}
