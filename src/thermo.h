#ifndef THERMO_H
#define THERMO_H

#include <fstream>
#include <string>

#include "particles.h"
#include "box.h"

class Thermo
{
/* Thermo class, with the write() method. It writes timestep,
   temperature, kinetic, potential and total energy */

  std::ofstream *thermofile;
  void write_header();

 public:
  int nfreq; /* How often we write the thermo information*/
  Thermo(int nfreq, std::ofstream *fname);
  ~Thermo();
  void write(int step, Particles *part);
};

#endif
