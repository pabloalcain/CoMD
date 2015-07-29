#include "cell.h"

Cell::Cell(int nidx) {
  natoms = 0;
  /* cell density < 2x avg. density */
  idxlist = (int *) malloc(nidx*sizeof(int));
};

void Cell::add_atom(int idx) {
  /* Add the atom idx in position "natoms" */
  idxlist[natoms] = idx;
  natoms++;
};

CellList::CellList(double length, Particles *part, Potential *pot, Box *box){
  /* length is the (max) size of the cell. It will be slightly changed
     to fit an integer number of cells in the simulation box. N is the
     total number of particles.

     I don't like that much passing N, since it is only used to find
     the maximum cell list*/
  
  /* Create the cell list */
  ncells = 1;
  for (int i=0; i<3; i++) {
    ngrid[i]  = floor(box->size[i] / length);
    delta[i]  = box->size[i] / ngrid[i];
    boxoffs[i] = 0.5 * box->size[i] - 0.5 * delta[i];
    ncells *= ngrid[i];
  }
  
  /* allocate cell list storage */
  list = (Cell *) malloc(ncells*sizeof(Cell));
  pairs = (int *) malloc(2*ncells*ncells*sizeof(int));

  /* allocate index lists within cell. cell density ~< 2x avg. density */
  nidx = 2*part->N / ncells + 2;
  for (int i=0; i<ncells; ++i) {
    list[i] = Cell(nidx);
  }

  /* build cell pair list, assuming newtons 3rd law. */
  npairs = 0;
  for (int i=0; i < ncells-1; ++i) {
    int i1[3];

    i1[0] = i % ngrid[0];
    i1[1] = i % (ngrid[0] * ngrid[1]) / ngrid[0];
    i1[2] = i / (ngrid[0] * ngrid[1]);

    for (int j=i+1; j<ncells; ++j) {
      int i2[3];
      double r[3];
      
      i2[0] = j % ngrid[0];
      i2[1] = j % (ngrid[0] * ngrid[1]) / ngrid[0];
      i2[2] = j / (ngrid[0] * ngrid[1]);

      for (int k=0; k < 3; k++) {
	r[k] = (i2[k] - i1[k]) * delta[k];
      }

      /* check for cells on a line that are too far apart */
      for (int k=0; k < 3; k++) {
	if (fabs(r[k]) > pot->rcut + delta[k]) continue;
      }
      /* check for cells in a plane that are too far apart */

      for (int k=0; k < 3; k++) {
	int k2 = (k + 1) % 3;
	double dist = r[k]*r[k] + r[k2]*r[k2];	
	double delta2 = delta[k]*delta[k] + delta[k2]*delta[k2];
	if (sqrt(dist) > pot->rcut + sqrt(delta2)) continue;
      }

      /* other cells that are too far apart */
      for (int k=0; k<3; k++) {
	double dist = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	double delta3 = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
	if (sqrt(dist) > pot->rcut + sqrt(delta3)) continue;
      }
      /* cells are close enough. add to list */
      pairs[2*npairs+0] = i;
      pairs[2*npairs+1] = j;
      npairs++;
    }
  }
  std::cout << "Cell list has " << ngrid[0] << "x" << ngrid[1] << "x" << ngrid[2]
	    << "=" << ncells << " cells with" << npairs << "/" << ncells*(ncells-1)/2
	    << " pairs and " << nidx << " atoms/celllist." << std::endl;
};

void CellList::update(Particles *part)
{
  for (int i=0; i < ncells; i++) {
    list[i].natoms=0;
  }

  int maxidx=0;
  for (int atom=0; atom < part->N; atom++) {
    int i1[3];
    int cellidx;
    for (int j= 0; j < 3; j++)
      i1[j] = part->x[atom*3 + j]/delta[j];

    cellidx = ngrid[0]*ngrid[1]*i1[2]+ngrid[0]*i1[1]+i1[0];

    list[cellidx].add_atom(atom);
    int idx = list[cellidx].natoms;
    if (idx > maxidx) maxidx = idx;
  }
  if (maxidx > nidx) {
    std::cout << "Error: overflow in cell list: " << maxidx << "/" << nidx
	      << " atoms/cell." << std::endl;
    exit(1);
  }
  return;
}
