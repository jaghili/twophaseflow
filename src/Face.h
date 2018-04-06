#ifndef FACE_H
#define FACE_H

#include "defs.h"

struct Face {
  uint idx;
  
  uint iT_left; // indice element gauche
  uint iT_right; // indice element droite
  
  Vec xa;  // only for plotting purposes
  Vec xb;
  Vec xF; // centre de la face
  R vol; // length of the face

  R flux_oil=0.; // flux
  R flux_wat=0.;

  Vec nF;
};

#endif
