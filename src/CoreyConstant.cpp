#include "CoreyConstant.h"

CoreyConstant::CoreyConstant() { 
  vg = Vec::Zero(dim);
  vg(0) = gx;
  vg(1) = gy;

  // Capillary pressures  
  pcname="corey-constant";
  _pc1 = 0.; // do not touch !
  _pc2 = (rhow-rhoo) * grav * L;
  b1 = 0.; // b1 < b2
  b2 = _pc2 ; // assert pc2 >= b2
  tau1 = 1. - b1/_pc2;
  tau2 = tau1 + 1 + b1 * std::log(1-tau1)/_pc2;
  tau3 = tau2 + 1;

  
  // taumin
  taumin = 0.;
  taumax = tau3-1e-12;
      
  // capillary pressures
  pc[0][0] = [this](R sat) -> R { return _pc1; }; // constant dans frac
  pc[1][1] = [this](R sat) -> R {
    return _pc2 - b2 * std::log(1-sat); // Corey dans matrice
  };
  pc[0][1] = [this](R tau) -> R {
    R s;
    if (tau < tau1) { s = _pc1; }
    else if (tau <= tau2) {
      s = _pc2*(tau-tau1); 
    }
    else {
      s = _pc2 - b2 * std::log(1-(tau-tau2));
    }
    return s;
  }; // pc01
  pc[1][0] = pc[0][1];
    
  // derivative
  dpc[0][0] = [this](R sat) -> R { return 0.; }; // pc00
  dpc[1][1] = [this](R sat) -> R { return b2/(1.-sat); }; // pc11
  dpc[0][1] = [this](R tau) -> R {
	R s;
	if (tau < tau1) { s = 0.; }
	else if (tau <= tau2) { s = _pc2;}
	else { s =  b2/(1-(tau-tau2)); }
	return s;
  }; // pc01
  dpc[1][0] = dpc[0][1];

  // saturations
  sat[0][0] = [this](R tau) -> R { return 0.; };
  sat[1][1] = [this](R tau) -> R { return 0.; };
  sat[0][1] = [this](R tau) -> R {
    R s;
    if (tau<= tau1) {
      s = tau;
    }
    else {
      s = 1;
    }
    return s;
  };
  sat[1][0] = [this](R tau) -> R {
    return (tau>=tau2)*(tau-tau2);
  };

  
  dsat[0][0] = [this](R tau) -> R { return 0.; };
  dsat[0][1] = [this](R tau) -> R {
    R s;
    if (tau<= tau1) { s = 1.; }
    else { s =  0.;  }
    return s;
  };
  dsat[1][0] = [this](R tau) -> R { return (tau>tau2)*1.; };
  dsat[1][1] = [this](R tau) -> R { return 0.; };
};
