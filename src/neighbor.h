#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include "particles.h"
#include "cell.h"
#include "potential.h"

class Neighbor
{
  friend class System;
  friend class CellList;
  friend class CoMD; 
 int *list;
  int *num;
  double **cutoff, maxcutoff;
  int nfreq;

  public:
  Neighbor(Particles *part, Potential *pot, double skin, int nfreq);
  void update(CellList *cells, Particles *part, Potential *pot, Box *box);
};

#endif
