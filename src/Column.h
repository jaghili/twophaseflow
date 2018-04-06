#ifndef COLUMN_HPP
#define COLUMN_HPP

#include "Mesh.h"
#include "Physics.h"

struct Column : Mesh {
  Column(int n, double w, Physics & P); // Build a column 
};

#endif
