#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "Mesh.h"
#include "Data.h"
#include "Physics.h"
#include "Residual.h"

// Assemble the real jacobian
SpMat computeJacobian(Mesh* Mh, Dat & X, Dat & Xold, Physics & P);

// Approx. jacobian: requires flux computation (slow)
Mat computeApproxJacobian (Mesh* Mh, Dat & X, Dat & Xold, Physics & P);
#endif
