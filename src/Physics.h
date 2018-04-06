#ifndef PHYSICS_H
#define PHYSICS_H

#include "Data.h"

struct Physics {
  const uint dim = 2;
  
  // domain bounds
  const float a = 0;
  const float b = 100.;

  const float rhow = 1e3;
  const float rhoo = 800;

  const float L = b-a;

  float grav = 10;
  const float gx = 0;
  const float gy = -grav;

  Vec vg = Vec::Zero(dim);
  
  const float sor = 0.0;
  const float swi = 0.0;
  const float muo = 5e-3;
  const float muw = 1e-3;

  const R pw = 0.; // base pressure water 
  
  const R phi = 0.1; // porosity

  static const short Nrt = 2;
  
  const R perm[Nrt] = {
    1e-13, // perm of rt1
    1e-15 // perm of rt2
  };
  
  // boundary
  R s_bottom = 1.; // oil
  const R s_top = 0.; // oil
  const R pw_inj = 1e5; // water injection 
  const R pw_bottom = pw_inj + L*grav*rhow; // water bottom
  const R pw_top = pw;  // water top

  // Mobilites
  R Mo(R so) { return (so>=sor) * std::pow((so-sor)/(1.-sor),2) / muo; }
  R Mw(R sw) { return (sw>=swi) * std::pow((sw-swi)/(1.-swi),2) / muw; }
  R dMo(R so) { return (so>=sor) * 2 * (so-sor)/std::pow((1-sor),2) / muo; }
  R dMw(R sw) { return (sw>=swi) * 2 * (sw-swi)/std::pow((1-swi),2) / muw; }

  // R Mo(R so) { return (so>=sor) * (so-sor)/(1.-sor) / muo; }
  // R Mw(R sw) { return (sw>=swi) * (sw-swi)/(1.-swi) / muw; }
  // R dMo(R so) { return (so>=sor) /(1-sor) / muo; }
  // R dMw(R sw) { return (sw>=swi) /(1-swi) / muw; }


  // project données
  float eps1=1e-6;
  
  void project(Dat & X, const Dat & Xold, const R & tau1, const R & tau2, const R & tau3) {
    // project saturations
    for (int i = 0; i < X.S.rows(); i++) {
      X.S(i) = (X.S(i)<0 ? 0. : X.S(i));
      X.S(i) = (X.S(i)>1. ? 1. - 1e-12 : X.S(i));
    }
    // project tau's
    for (int i = 0; i < X.Tau.rows(); i++) {
      if (X.Tau(i) > tau1 and Xold.Tau(i) < tau1 - eps1) {
	X.Tau(i) = tau1 - 0.5*eps1;
      }
      if (X.Tau(i) < tau1 and Xold.Tau(i) > tau1 + eps1) {
	X.Tau(i) = tau1 + 0.5*eps1;
      }
      if (X.Tau(i) > tau2 and Xold.Tau(i) < tau2 - eps1) {
	X.Tau(i) = tau2 - 0.5*eps1;
      }
      if (X.Tau(i) < tau2 and Xold.Tau(i) > tau2 + eps1) {
	X.Tau(i) = tau2 + 0.5*eps1;
      }
      if (X.Tau(i) < 0.) {
	X.Tau(i) = 1e-12;
      }
      if (X.Tau(i) >= tau3) {
	X.Tau(i) = tau3 - 1e-12;
      }
    }
  }
  
  // Capillary pressures  
  std::string pcname="pc";
  R _pc1; // do not touch !
  R _pc2;
  R b1; // b1 < b2
  R b2; // assert pc2 >= b2
  R tau1;
  R tau2;
  R tau3;

  // taumin
  R taumin;
  R taumax;
  
  // laws  
  std::function<R(R)> pc[Nrt][Nrt];
  std::function<R(R)> dpc[Nrt][Nrt];
  std::function<R(R)> sat[Nrt][Nrt];
  std::function<R(R)> dsat[Nrt][Nrt];

  // IO
  void writetofile();  
}; // end Physics

#endif
