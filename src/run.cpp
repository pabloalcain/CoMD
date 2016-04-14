#include "run.h"

int main(int argc, char* argv[]) {
  int nsteps;
  Box *box = new Box(20.0);
  Particles *particles = new Particles(20, 20, box);
  //Particles *particles = new Particles("init.lammpstrj", box);
  Potential *potential = new Potential(particles);
  //Integrator *integrator = new Integrator(0.1);
  Integrator *integrator = new Cooldown(0.1, 0.9, 200);
  
  
  std::ofstream file;
  file.open("dump.lammpstrj");
  Dump *dump = new Dump(1, &file);
  
  std::ofstream file2;
  file2.open("thermo.out");
  Thermo *thermo = new Thermo(1, &file2);
  
  if (argc != 2) {
    std::cerr << "usage: " << argv[0] << " nsteps" << std::endl;
    exit(1);
  }
  
  nsteps = atoi(argv[1]);
  std::cout << "Initializing system..." << std::endl;
  System sys = System(box, particles, potential, integrator, dump, thermo);
  sys.run(nsteps);
}
