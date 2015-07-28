#include "comd.h"

int main(int argc, char* argv[]) {
    Box box = Box(5.0);
    Particles particles = Particles(500, 500, box);
    Potential potential = Potential();
    Integrator integrator = Integrator(0.01);

    std::cout << "Initializing system" << std::endl;
    System sys = System(&box, &particles, &potential, &integrator);
    sys.run(1000);
}
