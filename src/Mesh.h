#ifndef MESH_H
#define MESH_H

#include "defs.h"
#include "Element.h"
#include "Face.h"
#include "BFace.h"
#include "IFace.h"

struct Mesh {
  uint d;

  uint N; // nombre de mailles
  uint Ni; // nombre d'interfaces avec diff rt
  uint Nf; // nombre de faces standards
  uint Nfb; // nombre de faces de bord

  uint Nrt;

  double xmin=0.;
  double xmax=0.;
  double ymin=1.;
  double ymax=1.;
  
  std::vector<Element> Mailles; // toutes les mailles

  std::vector<Face> FacesStd; // faces standards (same rt)
  std::vector<IFace> Interfaces; // faces interfaces (diff rt)
  std::vector<BFace> BoundaryFaces; // faces de bords

  bool ELEMENTS_ARE_SORTED = false;
  
  //Mesh(); // default ctor
  void setRocktype(uint _idx, uint rocktype); // set pointwise rt
  void computeInterfaces(); // once rocktype is fixed, compute interfaces and fill intefaces vector
  void display();
};

#endif
