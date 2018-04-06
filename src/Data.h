#ifndef DATA_H
#define DATA_H

#include "defs.h"
#include "Mesh.h"

struct Dat {
  uint N; // get from Mesh
  uint Ni; // get from Mesh
  uint Ndof; // N+Ni
  
  Vec S;
  Vec P;
  Vec Pi;
  Vec Tau;
  Dat(Mesh* Mh);
  Vec get();
  void set(Vec & U);
  Vec gradP();
  void display();
};
#endif
