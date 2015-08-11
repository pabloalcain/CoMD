#include "run.h"

int main(int argc, char* argv[]) {
    int nsteps;
    Box box = Box(20.0);
    Particles particles = Particles(20, 20, box);
    Potential potential = Potential();
    Integrator integrator = Integrator(0.0005, 100);

    std::ofstream file;
    file.open("dump.lammpstrj");
    Dump dump = Dump(100, &file);

    std::ofstream file2;
    file2.open("thermo.out");
    Thermo thermo = Thermo(1, &file2);
    
    if (argc != 2) {
	std::cerr << "usage: " << argv[0] << " nsteps" << std::endl;
	exit(1);
    }
    nsteps = atoi(argv[1]);
    std::cout << "Initializing system..." << std::endl;
    System sys = System(&box, &particles, &potential, &integrator, &dump, &thermo);
    sys.run(nsteps);
}
