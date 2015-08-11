#include "run.h"

int main(int argc, char* argv[]) {
    int nsteps;
    Box *box = new Box(20.0);
    Particles *particles = new Particles(20, 20, *box);
    Potential *potential = new Potential();
    Integrator *integrator = new Cooldown(0.001, 0.99, 100);

    std::ofstream file;
    file.open("dump.lammpstrj");
    Dump *dump = new Dump(100, &file);

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
