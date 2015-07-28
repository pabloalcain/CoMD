#ifndef BOX_H
#define BOX_H

class Box
{
  /* Box structure. Needed to create particles in it */
  friend class Particles;
  double size[3];
  
 public:
  Box(double L);
};
#endif

