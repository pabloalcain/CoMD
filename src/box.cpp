#include "box.h"

Box::Box(double L) {
  for (int i=0; i<3; i++)
    size[i] = L;
}
