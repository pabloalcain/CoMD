#include "dump.h"

Dump::Dump(int _nfreq, std::ofstream *fname) {
  dumpfile = fname;
  nfreq = _nfreq;
}

Dump::~Dump() {
  dumpfile->close();
}

void Dump::write(int step, Particles *part, Box *box) {
  write_header(step, part, box);
  write_part(part);
}

void Dump::write_header(int step, Particles *part, Box *box) {
  *dumpfile << "ITEM: TIMESTEP" << std::endl;
  *dumpfile << step << std::endl;
  *dumpfile << "ITEM: NUMBER OF ATOMS" << std::endl;
  *dumpfile << part->N << std::endl;
  *dumpfile << "ITEM: BOX BOUNDS p p p" << std::endl;
  *dumpfile << -box->size[0]/2 << " " << box->size[0]/2 << std::endl;
  *dumpfile << -box->size[1]/2 << " " << box->size[1]/2 << std::endl;
  *dumpfile << -box->size[2]/2 << " " << box->size[2]/2 << std::endl;
  *dumpfile << "ITEM: ATOMS id type x y z vx vy vz occ" << std::endl;
}

void Dump::write_part(Particles *part) {
  for (int i = 0; i < part->N; i++) {
    int is = part->isospin[i]?2:1;
    *dumpfile << i;
    *dumpfile << " " << is;
    for (int l = 0; l < 3; l++) 
      *dumpfile << " " << part->x[3*i + l]; 
    for (int l = 0; l < 3; l++) 
      *dumpfile << " " << part->v[3*i + l];
    *dumpfile << " " << part->occ[i];
    *dumpfile << std::endl;
  }
}
