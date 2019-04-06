#include "Constants.h"

///////////////////////////////////////
// LOIS "CONSTANTES" 
///////////////////////////////////////

Constants::Constants() { 
  vg = Vec::Zero(dim);
  vg(0) = gx;
  vg(1) = gy;

  // Capillary pressures  
  pcname="constant";
  pef = 0.; // do not touch !
  pem = (rhow-rhoo) * grav * L / 10.;
  b1 = 0.; // b1 < b2
  b2 = 0.; // assert pc2 >= b2
  tau1 = 1;
  tau2 = 2;
  tau3 = 3;

  // range of tau
  taumin = 0.;
  taumax = tau3-1e-12;
      
  // capillary pressures
  pc[0][0] = [this](R sat) -> R { return pef; };
  pc[1][1] = [this](R sat) -> R { return pem; };
  pc[0][1] = [this](R tau) -> R {
    return (tau<=1)*pef
    +(1.<tau and tau<=2.) * (pef +(tau-tau1)*(pem-pef))
    + (tau>2.)*pem;
  }; // pc01
  pc[1][0] = pc[0][1];

    
  // derivative
  dpc[0][0] = [this](R sat) -> R { return 0.; }; // pc00
  dpc[1][1] = [this](R sat) -> R { return 0.; }; // pc11
  dpc[0][1] = [this](R tau) -> R {
    return (tau<=tau2)*0. + (1<tau && tau<=2)*(pem-pef) + (tau>2.)*0.;
  }; // pc01
  dpc[1][0] = dpc[0][1];

  // saturations
  sat[0][0] = [this](R tau) -> R { return 0.; };
  sat[1][1] = [this](R tau) -> R { return 0.; };
  sat[0][1] = [this](R tau) -> R {
    return (tau<=1)*tau + (1<tau && tau<=2)*1. + (tau>2)*1.;
  };
  sat[1][0] = [this](R tau) -> R {
    return (tau<=1)*0. + (1<tau && tau<=2)*0. + (tau>2)*(tau-2.);
  };

  
  dsat[0][0] = [this](R tau) -> R { return 0.; };
  dsat[0][1] = [this](R tau) -> R {
    return (tau<=1)*1. + (1<tau && tau<=2)*0. + (tau>2)*0.; 
  };
  dsat[1][0] = [this](R tau) -> R { return (tau<=1)*0. + (1<tau && tau<=2)*0. + (tau>2)*1.; };
  dsat[1][1] = [this](R tau) -> R { return 0.; };
};
