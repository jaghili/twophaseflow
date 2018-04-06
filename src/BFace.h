#ifndef BFACE_H
#define BFACE_H

#include "Face.h"

struct BFace : Face {
  bool DIRICHLET_BC = false;
  R sat; // Dirichlet saturation
  R pw; // Dirichlet pressure

  bool NEUMANN_BC = false;
  R neumann_flux_oil; // Neumann condition
  R neumann_flux_wat; // Neumann condition
  
  BFace();
  void setNeumannBC(R _so, R _pw, R _fo, R _fw); // Neumann
  void setDirichletBC(R _so, R _pw); // Dirichlet
};

#endif
