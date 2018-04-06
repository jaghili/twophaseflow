#ifndef FVCA5MESH_H
#define FVCA5MESH_H

#include "Mesh.h"

struct FVCA5Mesh : Mesh {
  FVCA5Mesh(const char* file, R scale);
};

#endif
