#include "neighbor.h"

Neighbor::Neighbor(Particles *part, Potential *pot, double skin, int _nfreq) {
  list = (int *) malloc(part->N * part->N * sizeof(int));
  num = (int *) malloc(part->N * sizeof(int));

  maxcutoff = -1.0;
  int n = part->ntypes;
  cutoff = (double **) malloc((n+1)*sizeof(double *));
  for (int i = 0; i < n+1; i++) {
    cutoff[i] = (double *) malloc((n+1)*sizeof(double));
    for (int j = 0; j < n+1; j++) {
      cutoff[i][j] = pot->rcut[i][j] > 3 * part->sigma_r? pot->rcut[i][j] : 3 * part->sigma_r;
      cutoff[i][j] += skin;
      if (cutoff[i][j] > maxcutoff) maxcutoff = cutoff[i][j];
    }
  }
  nfreq = _nfreq;
}

void Neighbor::update(CellList *cells, Particles *part, Potential *pot, Box *box) {
  /* Add the atom idx in position "natoms" */
  double x[3];
  double dr;
  int n;
  int t1, t2;

  cells->update(part, box);
  
  for (int i = 0; i < part->N; i++) num[i] = 0;
  for (int k = 0; k < cells->ncells; k++) {
    Cell *cell = cells->list + k;
    for (int i = 0; i < cell->natoms - 1; i++) {
      int ii = cell->idxlist[i];
      n = 0;
      x[0] = part->x[3*ii+0];
      x[1] = part->x[3*ii+1];
      x[2] = part->x[3*ii+2];
      t1 = part->isospin[ii]?2:1;
      for (int j = i+1; j < cell->natoms; j++) {
        t2 = part->isospin[ii]?2:1;
        int jj = cell->idxlist[j];
        dr = 0;
        for (int l = 0; l < 3; l++)
          dr += box->pbc(x[l] - part->x[3*jj + l], l) * box->pbc(x[l] - part->x[3*jj + l], l);
        if (dr < cutoff[t1][t2] * cutoff[t1][t2]){
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
      t1 = part->isospin[ii]?2:1;
      n = num[ii];
      x[0] = part->x[3*ii+0];
      x[1] = part->x[3*ii+1];
      x[2] = part->x[3*ii+2];
      for (int j = 0; j < c2->natoms; j++) {
        int jj = c2->idxlist[j];
        t2 = part->isospin[ii]?2:1;
        dr = 0;
        for (int l = 0; l < 3; l++)
          dr += box->pbc(x[l] - part->x[3*jj + l], l) * box->pbc(x[l] - part->x[3*jj + l], l);

        if (dr < pot->rcut[t1][t2] * pot->rcut[t1][t2]){
          list[ii*part->N + n] = jj;
          n++;
        }
      }
      num[ii] = n;
    }
  }
}
