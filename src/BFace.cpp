#include "BFace.h"

BFace::BFace() {
  flux_oil = 0.;
  flux_wat = 0.;
  pw = 0.;
  sat = 0.;
}

void BFace::setNeumannBC(R _so, R _pw, R _fo, R _fw) { // Neumann
  DIRICHLET_BC = false;
  NEUMANN_BC = true;
  neumann_flux_oil = _fo;
  neumann_flux_wat = _fw;
  sat = _so;
  pw = _pw;
}

void BFace::setDirichletBC(R _so, R _pw) { // Neumann
  DIRICHLET_BC = true;
  NEUMANN_BC = false;
  sat = _so;
  pw = _pw;
}
