#include "Corey.h"

///////////////////////////////////////
// LOIS DE COREY TESTEES ET VALIDEES //
///////////////////////////////////////

Corey::Corey() { 
  vg = Vec::Zero(dim);
  vg(0) = gx;
  vg(1) = gy;

  // Capillary pressures  
  pcname="corey";
  pef = 0.;
  pem = (rhow-rhoo) * grav * L / 2;
  b1 = pem/10.; // b1 < b2
  b2 = pem; // assert pc2 >= b2
  tau1 = 1.0 - b1/pem;
  tau2 = tau1 + 1 + b1 * std::log(1-tau1)/pem;
  tau3 = tau2 + 1;

  // taumin
  taumin = 0.;
  taumax = tau3-1e-12;
  
  // capillary pressures
  pc[0][0] = [this](R sat) -> R { return pef - b1 * std::log(1.-sat); }; // frac rt0
  pc[1][1] = [this](R sat) -> R { return pem - b2 * std::log(1.-sat); }; // mat rt1
  pc[0][1] = [this](R tau) -> R {
    R s;
    if (tau < tau1) { s = pef - b1 * std::log(1-tau); }
    else if (tau <= tau2) {
      s = pef - b1*std::log(1-tau1) + (tau-tau1)*(pem-pef);
    }
    else { s = pem - b2 * std::log(1-(tau-tau2)); }
    return s;
  }; // pc01
  pc[1][0] = pc[0][1];
    
  // derivative
  dpc[0][0] = [this](R sat) -> R { return b1/(1.-sat); }; // pc00
  dpc[1][1] = [this](R sat) -> R { return b2/(1.-sat); }; // pc11
  dpc[0][1] = [this](R tau) -> R {
	R s;
	if (tau < tau1) { s = b1 / (1-tau); }
	else if (tau <= tau2) { s = pem-pef;}
	else { s =  b2 / (1-tau+tau2); }
	return s;
  }; // pc01
  dpc[1][0] = dpc[0][1];

  // saturations
  sat[0][0] = [this](R tau) -> R { return 0.; };
  sat[1][1] = [this](R tau) -> R { return 0.; };
  sat[0][1] = [this](R tau) -> R {
    R s;
    if (tau <= tau1) { s = tau; }
    else if ( tau <= tau2) {
      s = 1 - (1-tau1) * std::exp(-(tau-tau1) * pem/b1);
    }
    else { s = 1 - std::pow(1-(tau-tau2),b2/b1) * std::exp(-pem/b1); }
    return s;
  };
  sat[1][0] = [this](R tau) -> R {
    return (tau>tau2)*(tau-tau2);
  };

  
  dsat[0][0] = [this](R tau) -> R { return 0.; };
  dsat[0][1] = [this](R tau) -> R {
    R s;
    if (tau<= tau1) {       s = 1.;    }
    else if (tau <= tau2) { s = pem/b1 * (1-tau1) * std::exp(-(tau-tau1)*pem/b1); }
    else { s =  (b2/b1) * std::pow(1-(tau-tau2),b2/b1 - 1)*std::exp(-pem/b1);  }
    return s;
  };
  dsat[1][0] = [this](R tau) -> R { return (tau>tau2)*1.; };
  dsat[1][1] = [this](R tau) -> R { return 0.; };
};
