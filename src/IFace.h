#ifndef IFACE_H
#define IFACE_H

#include "Face.h"

struct IFace : Face {
  // volume fictiv
  uint v_iT_left;
  uint v_iT_right;
  
  // fluxes
  R flux_oil_left=0.;
  R flux_oil_right=0.;

  R flux_wat_left=0.;
  R flux_wat_right=0.;

  bool is_taumin = false; // default tau
  R tau = -1e300;
};

#endif
