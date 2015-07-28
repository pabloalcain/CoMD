#ifndef CELL_H
#define CELL_H

#include "sys.h"
class System;
class Cell
{
  /* Cell structure */
  friend class Particles;
  friend class Box;

  int natoms; /* number of atoms in this cell */
  int *idxlist; /* list of atom indices in the cell */
  
 public:
  Cell(int idx);
  void add_atom(int idx);
};

class CellList
{
  /* List of cells */
  friend class Cell;
  int nidx; /* Maximum number of atoms indices in each cell */
  int ngrid[3]; /* Number of cells in each dimension */
  int delta[3]; /* Actual size of the cells in each dimension */
  int boxoffs[3]; /* Box offset in each dimension */

  int ncells; /* Total number of cells */
  Cell *list; /* Actual list of cells */

  int npairs; /* Number of pairs that interact */
  int *pairs; /* Pairs of cells that are within cutoff distance */
  
 public:
  CellList(double length);
  void update(System *sys);
};
   
#endif
