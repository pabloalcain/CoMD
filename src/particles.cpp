#include "particles.h"

Particles::Particles(int nprot, int nneut, Box box){
    N = nprot + nneut;
    x = (double *) malloc(3 * N * sizeof(double));
    v = (double *) malloc(3 * N * sizeof(double));
    f = (double *) malloc(3 * N * sizeof(double));
    g = (double *) malloc(3 * N * sizeof(double));
    isospin = (bool *) malloc(N * sizeof(char));
    spin = (bool *) malloc(N * sizeof(char));
    
    int L = (int) ceil(pow(N, 1.0/3.0));
    double dx = box.size[0]/L;
    double dy = box.size[1]/L;
    double dz = box.size[2]/L;
    
    int idx;
    /* By default we put the particles on a lattice */
    
    for (int ix = 0; ix < L; ix++) {
	for (int iy = 0; iy < L; iy++) {
	    for (int iz = 0; iz < L; iz++) {
		idx = ix * L * L + iy * L + iz;
		if (idx >= N) continue;
		x[3*idx + 0] = ix * dx - box.size[0]/2;
		x[3*idx + 1] = iy * dy - box.size[1]/2;
		x[3*idx + 2] = iz * dz - box.size[2]/2;
	    }
	}
    }

    for (int i = 0; i < 3*N; i++){
	v[i] = 0.5*(double)rand()/RAND_MAX;
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
    sig_r = 1.0;
    sig_p = 1.0;
    mass = 1.0;
}
