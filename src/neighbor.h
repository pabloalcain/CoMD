#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include "particles.h"
#include "cell.h"

class Neighbor
{
  friend class System;
  int *list;
  int *num;
  double skin;

  public:
  Neighbor(Particles *part);
  void update(CellList *cells, Particles *part, Potential *pot, Box *box);
};

#endif
