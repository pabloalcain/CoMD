#ifndef DUMP_H
#define DUMP_H

#include <fstream>
#include <string>

#include "particles.h"
#include "comd.h"
#include "box.h"

class Dump
{
/* Dump class, with the write() method. Prepared to mimic lammps
   style with id, type, x, y, z, vx, vy, vz */

  std::ofstream *dumpfile;
  void write_header(int step, Particles *part, Box *box);
  void write_part(Particles *part, CoMD *comd);
  
 public:
  int nfreq; /* How often we write the dump */
  Dump(int nfreq, std::ofstream *fname);
  ~Dump();
  void write(int step, Particles *part, CoMD *comd, Box *box);
};

#endif
