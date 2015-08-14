#include "particles.h"

void Particles::allocate() {
  x = (double *) malloc(3 * N * sizeof(double));
  v = (double *) malloc(3 * N * sizeof(double));
  f = (double *) malloc(3 * N * sizeof(double));
  occ = (double *) malloc(N * sizeof(double));
  isospin = (bool *) malloc(N * sizeof(char));
  spin = (bool *) malloc(N * sizeof(char));
}


Particles::Particles(int nprot, int nneut, Box *box){
  /* isospin false = type 1 = neutron
     isospin true = type 2 = proton */
  N = nprot + nneut;
  ntypes = 2;
  allocate();
   
  int L = cbrt(N);
  if (N != L * L * L) L++;
  double dx = box->size[0]/L;
  double dy = box->size[1]/L;
  double dz = box->size[2]/L;
    
  int idx;
  /* By default we put the particles on a lattice */
    
  for (int ix = 0; ix < L; ix++) {
    for (int iy = 0; iy < L; iy++) {
      for (int iz = 0; iz < L; iz++) {
	idx = ix * L * L + iy * L + iz;
	if (idx >= N) continue;
	x[3*idx + 0] = ix * dx - box->size[0]/2;
	x[3*idx + 1] = iy * dy - box->size[1]/2;
	x[3*idx + 2] = iz * dz - box->size[2]/2;
      }
    }
  }
  double vcm[3];
  for (int j = 0; j < 3; j++) vcm[j] = 0;
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < 3; j++) {
      v[3*i+j] = 0.04*(double)rand()/RAND_MAX;
      vcm[j] += v[3*i+j]/N;
    }
  }
  
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < 3; j++) {
      v[3*i+j]-= vcm[j];
    }
  }

  for (int i = 0; i < N; i++)
    occ[i] = 0.0;

    
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
  sigma_r = 1.3;
  sigma_p = 0.47;
  mass = 938.0;
}

Particles::Particles(const std::string& fname, Box *box){
  std::string line;
  std::ifstream dumpfile(fname.c_str());

  double number;
  int nline = 0;
  int maxlines = -1;
  ntypes = 2;

  while(getline(dumpfile, line)) {
    nline++;
    if (nline == 4) {
      N = atoi(line.c_str());
      maxlines = N + 9;
      allocate();
    }
    if (nline >= 6 && nline < 9) {
      std::istringstream iss(line);
      iss >> number;
      double size = - number;
      iss >> number;
      size += number;
      box->size[nline - 6] = size;
    }
    if (nline >= 10) {
      std::istringstream iss(line);
      int idx;
      int type;
      iss >> idx;
      iss >> type;
      isospin[idx] = (type == 2);
      iss >> type;
      spin[idx] = (type == 1);
      for (int i = 0; i < 3; i++) {
	iss >> number;
	x[3*idx + i] = number;
      }
      for (int i = 0; i < 3; i++) {
	iss >> number;
	v[3*idx + i] = number;
      }
    }
    if (nline == maxlines) break;
  }
  /* Dummy values */
  sigma_r = 1.3;
  sigma_p = 0.47;
  mass = 938.0;
}
