#ifndef SYS_H
#define SYS_H

#include <math.h>
#include <iostream>


#include "particles.h"
#include "potential.h"
#include "integrator.h"
#include "cell.h"
#include "box.h"
#include "dump.h"

class System {
  /* This is the basic system structure */
  Box *box;
  Particles *part;
  Potential *pot;
  Integrator *integ;
  CellList *cells;
  Dump *dump;
 public:
  System(Box *_box, Particles *_part, Potential *_pot, Integrator *_integ, Dump *_dump);
  void forces();
  void forces_all();
  void run(int nsteps);
};

#endif
