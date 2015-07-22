#include "particles.h"

Particles::Particles(int nprot, int nneut, Box box){
    N = nprot + nneut;
    x = (double *) malloc(3 * N * sizeof(double));
    v = (double *) malloc(3 * N * sizeof(double));
    f = (double *) malloc(3 * N * sizeof(double));
    isospin = (bool *) malloc(N * sizeof(char));
    spin = (bool *) malloc(N * sizeof(char));
    
    int L = (int) ceil(pow(N, 1.0/3.0));
    double dx = box.Lx/L;
    double dy = box.Ly/L;
    double dz = box.Lz/L;
    
    int idx;
    /* By default we put the particles on a lattice */
    
    for (int ix = 0; ix < L; ix++) {
	for (int iy = 0; iy < L; iy++) {
	    for (int iz = 0; iz < L; iz++) {
		idx = ix * L * L + iy * L + iz;
		if (idx >= N) continue;
		x[3*idx + 0] = ix * dx;
		x[3*idx + 1] = iy * dy;
		x[3*idx + 2] = iz * dz;
	    }
	}
    }
    
    /* We set the first as neutrons, the last as protons */
    for (int i = 0; i < nneut; i++)
	isospin[i] = false;
    for (int i = nneut; i < N; i++)
	isospin[i] = true;

    /* Now we alternate the spin */
    
    bool status = true;
    for (int i = 0; i < N; i++){
	status = !status;
	spin[i] = status;
    }
 

    /* Dummy values */
    double sig_r = 1.0;
    double sig_p = 1.0;
    double mass = 1.0;
}
