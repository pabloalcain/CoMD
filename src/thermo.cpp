#include "thermo.h"

Thermo::Thermo(int _nfreq, std::ofstream *fname) {
  thermofile = fname;
  nfreq = _nfreq;
  write_header();
}

Thermo::~Thermo() {
  thermofile->close();
}

void Thermo::write(int step, Particles *part) {
  *thermofile << step << " " << part->ke*1.5 << " " << part->ke << " "
	      << part->pe << " " << part->ke + part->pe << std::endl;
}

void Thermo::write_header() {
  *thermofile << "Step Temperature Kinetic Potential Total" << std::endl;
}
