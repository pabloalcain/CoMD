#ifndef SYS_H
#define SYS_H

#include <math.h>
#include <iostream>


#include "particles.h"
#include "potential.h"
#include "integrator.h"
#include "cell.h"
#include "box.h"

class CellList;

class System {
  /* This is the basic system structure */
  Box *box;
  Particles *part;
  Potential *pot;
  Integrator *integ;
  CellList *cells;
 public:
  System(Box *_box, Particles *_part, Potential *_pot, Integrator *_integ);
  void forces();
  void run(int nsteps);
};

#endif
