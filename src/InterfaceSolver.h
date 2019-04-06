#ifndef INTERFACESOLVER_H
#define INTERFACESOLVER_H

#include "defs.h"
#include "Mesh.h"
#include "Data.h"
#include "Physics.h"
#include "BasicFunctions.h"
#include "IFace.h"

Vec InterfaceSolver(Mesh* Mh,
		    Dat & X,
		    IFace & iF,
		    Physics & P,
		    const R tolzero,
		    const R toldicho,
		    const R tolnwt
		    );

#endif
