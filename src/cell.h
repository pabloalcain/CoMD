#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <math.h>
#include <stdlib.h>

#include "box.h"
#include "particles.h"
#include "potential.h"


class Cell
{
  friend class CellList;
  friend class System;
  friend class CoMD;
  friend class Neighbor;
  /* Cell structure */
  int natoms; /* number of atoms in this cell */
  int *idxlist; /* list of atom indices in the cell */
  
 public:
  Cell(int nidx);
  void add_atom(int idx);
};

class CellList
{
  /* List of cells */
  friend class Cell;
  friend class System;
  friend class CoMD;
  friend class Neighbor;
  int nidx; /* Maximum number of atoms indices in each cell */
  int ngrid[3]; /* Number of cells in each dimension */
  int boxoffs[3]; /* Box offset in each dimension */
  double delta[3]; /* Actual size of the cells in each dimension */

  int ncells; /* Total number of cells */
  Cell *list; /* Actual list of cells */

  int npairs; /* Number of pairs that interact */
  int *pairs; /* Pairs of cells that are within cutoff distance */
  
 public:
  CellList(double length, Particles *part, double cutoff, Box *box);
  void update(Particles *part, Box *box);
};
   
#endif
