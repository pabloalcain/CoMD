#ifndef BOX_H
#define BOX_H

class Box
{
  /* Box structure. Needed to create particles in it */
  friend class Particles;
  friend class CellList;
  friend class Dump;
  friend class CoMD;
  double size[3];
  
 public:
  Box(double L);
  inline double pbc(double x, int i) {
    /* i is whether the direction is x(0), y(1) or z (2) */
    while (x >  size[i]/2) x -= size[i];
    while (x < -size[i]/2) x += size[i];
    return x;
  }
};
#endif

