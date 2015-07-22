#include "sys.h"

System::System(Box *_box, Particles *_part, Potential *_pot, Integrator *_integ){
    box = _box;
    part = _part;
    pot = _pot;
    integ = _integ;
}

void System::forces() {
    double x, y, z;
    double dx, dy, dz;
    double dr, dphi;

    for (int i = 0; i < part->N; i++) {
	x = part->x[3*i + 0];
	y = part->x[3*i + 1];
	z = part->x[3*i + 2];
 	for (int j = i; j < part->N; j++) {
	    dx = x - part->x[3*j + 0];
	    dy = y - part->x[3*j + 1];
	    dz = z - part->x[3*j + 2];
	    dr = dx * dx + dy * dy + dz * dz;
	    dr = sqrt(dr);
	    dphi = pot->dphi(dr);
	    part->f[3*i + 0] += dphi * dx;
	    part->f[3*i + 1] += dphi * dy;
	    part->f[3*i + 2] += dphi * dz;
	    /* Newton's 3rd law */
	    part->f[3*j + 0] -= dphi * dx;
	    part->f[3*j + 1] -= dphi * dy;
	    part->f[3*j + 2] -= dphi * dz;
	}
    }
}
    
void System::run(int nsteps) {
    for (int i = 0; i < nsteps; i++) {
	std::cout << "Step " << i << "out of "<< nsteps << std::endl;
	integ->first_step(part);
	forces();
	integ->final_step(part);
    }
}

