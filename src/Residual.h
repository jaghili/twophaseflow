#ifndef RESIDUAL_H
#define RESIDUAL_H

#include "defs.h"
#include "Mesh.h"
#include "Data.h"
#include "Physics.h"


// Compute Residual onces fluxes are computed 
Mat computeResidual(Mesh* Mh, Dat & X, Dat & Xold, Physics & P);

#endif
