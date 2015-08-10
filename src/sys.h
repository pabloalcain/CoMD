#ifndef SYS_H
#define SYS_H

#include <math.h>
#include <iostream>


#include "particles.h"
#include "potential.h"
#include "integrator.h"
#include "cell.h"
#include "comd.h"
#include "box.h"
#include "neighbor.h"
#include "dump.h"
#include "thermo.h"
#include "units.h"

class System {
  /* This is the basic system structure */
  Box *box;
  Particles *part;
  Potential *pot;
  Integrator *integ;
  CellList *cells;
  Dump *dump;
  Thermo *thermo;
  Units *units;
  CoMD *comd;
  Neighbor *neighbor;

 public:
  System(Box *_box, Particles *_part, Potential *_pot,
	 Integrator *_integ, Dump *_dump, Thermo *_thermo);
  void forces();
  void forces_all();
  void forces_neigh();
  void run(int nsteps);
};

#endif
