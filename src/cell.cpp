#include "cell.h"

Cell::Cell(int nidx) {
  natoms = 0;
  /* cell density < 2x avg. density */
  idxlist = (int *) malloc(nidx*sizeof(int));
}

void Cell::add_atom(int idx) {
  /* Add the atom idx in position "natoms" */
  idxlist[natoms] = idx;
  natoms++;
}

CellList::CellList(double length, Particles *part, double cutoff, Box *box){
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
  nidx = 10*part->N / ncells + 2;
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
      double dist, delta2, delta3;
      i2[0] = j % ngrid[0];
      i2[1] = j % (ngrid[0] * ngrid[1]) / ngrid[0];
      i2[2] = j / (ngrid[0] * ngrid[1]);


      for (int k=0; k < 3; k++) {
	r[k] = box->pbc((i2[k] - i1[k]) * delta[k], k);
      }

      /* check for cells on a line that are too far apart */
      if (fabs(r[0]) > cutoff + delta[0]) continue;
      if (fabs(r[1]) > cutoff + delta[1]) continue;
      if (fabs(r[2]) > cutoff + delta[2]) continue;

      /* check for cells in a plane that are too far apart */
      dist = r[0]*r[0] + r[1]*r[1];	
      delta2 = delta[0]*delta[0] + delta[1]*delta[1];
      if (sqrt(dist) > cutoff + sqrt(delta2)) continue;

      dist = r[2]*r[2] + r[1]*r[1];	
      delta2 = delta[2]*delta[2] + delta[1]*delta[1];
      if (sqrt(dist) > cutoff + sqrt(delta2)) continue;

      dist = r[0]*r[0] + r[2]*r[2];	
      delta2 = delta[0]*delta[0] + delta[2]*delta[2];
      if (sqrt(dist) > cutoff + sqrt(delta2)) continue;

      /* other cells that are too far apart */
      dist = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
      delta3 = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2];
      if (sqrt(dist) > cutoff + sqrt(delta3)) continue;

      /* cells are close enough. add to list */
      pairs[2*npairs+0] = i;
      pairs[2*npairs+1] = j;
      npairs++;
    }
  }
  std::cout << "Cell list has " << ngrid[0] << "x" << ngrid[1] << "x" << ngrid[2]
	    << "=" << ncells << " cells with " << npairs << "/" << ncells*(ncells-1)/2
	    << " pairs and " << nidx << " atoms/celllist." << std::endl;
}

void CellList::update(Particles *part, Box *box)
{
  for (int i=0; i < ncells; i++) {
    list[i].natoms=0;
  }
  
  int maxidx=0;
  for (int atom=0; atom < part->N; atom++) {
    
    int i1[3];
    double *pos = (part->x + atom*3);
    int cellidx;
    for (int j= 0; j < 3; j++){
      pos[j] = box->pbc(pos[j], j);
      i1[j] = floor((pos[j] + box->size[j]/2.0)/delta[j]);
    }
    cellidx = ngrid[0]*ngrid[1]*i1[2]+ngrid[0]*i1[1]+i1[0];
    if (cellidx > ncells - 1 || cellidx < 0) {
      std::cerr << "Error: cell " << i1[0] << ", " << i1[1] << ", " << i1[2] <<
	" does not exist. Maybe atom " << atom << ", x = " << pos[0] << ", " <<
	pos[1] << ", " << pos[2] << " is out of bounds?" << std::endl;
      exit(1);
    }
    
    list[cellidx].add_atom(atom);
    int idx = list[cellidx].natoms;
    if (idx > maxidx) maxidx = idx;
  }
  if (maxidx > nidx) {
    std::cerr << "Error: overflow in cell list: " << maxidx << "/" << nidx
	      << " atoms/cell." << std::endl;
    exit(1);
  }
  return;
}
