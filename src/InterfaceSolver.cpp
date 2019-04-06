#include "InterfaceSolver.h"
#include <iomanip>
#include <string>

// #define IS_DEBUG

extern uint nwt_max_iter;
extern uint nwt_is_fails;
extern uint nt;

// Interface solver algorithm
Vec InterfaceSolver(Mesh* Mh,
		    Dat & X,
		    IFace & F,
		    Physics & P,
		    const R tolzero,
		    const R toldicho,
		    const R tolnwt) {
  
  std::string prefix = "[IS" + std::to_string(F.idx) + "] ";

  // solution vector :  tau, ug, ul, pc
  Vec ResultVector = Vec::Zero(4); 
  
  R _tau = 0.;
  R _po =0.;
  R _pw =0.;
  
  // define left, right side
  uint left = F.iT_left; // K
  uint right = F.iT_right;  // L
  
  uint rtl = Mh->Mailles[left].rt;
  uint rtr = Mh->Mailles[right].rt;
  
  // ------------------------------------------------------------------------------------- //
  // oil
  R phio_left = X.P(left) + P.pc[rtl][rtl](X.S(left)) + P.rhoo * P.grav * Mh->Mailles[left].xT(1);
  R phio_right = X.P(right) + P.pc[rtr][rtr](X.S(right)) + P.rhoo * P.grav * Mh->Mailles[right].xT(1);

  // maille amont/aval
  uint o_amont = (phio_left >= phio_right ? left : right);
  uint o_aval = (phio_left >= phio_right ? right : left);

  // rocktypes
  uint rt_amont = Mh->Mailles[o_amont].rt;
  uint rt_aval = Mh->Mailles[o_aval].rt;
  assert(rt_amont != rt_aval);
  
  // phi amont, aval
  R phio_amont = X.P(o_amont) + P.pc[rt_amont][rt_amont](X.S(o_amont)) + P.rhoo * P.grav * Mh->Mailles[o_amont].xT(1);
  R phio_aval = X.P(o_aval) + P.pc[rt_aval][rt_aval](X.S(o_aval)) + P.rhoo * P.grav * Mh->Mailles[o_aval].xT(1);
  
  // Transmissibilities
  R Transo_amont = P.perm[rt_amont] / (F.xF - Mh->Mailles[o_amont].xT).norm();
  R Transo_aval = P.perm[rt_aval] / (F.xF - Mh->Mailles[o_aval].xT).norm();

  // Phio and its derivative
  std::function<R(R)> Phio = [&](R tau) -> R { // & capture all constants by ref
    if (P.Mo(X.S(o_amont)) > tolzero) {
      R u = P.Mo(X.S(o_amont)) * Transo_amont * phio_amont
	+ P.Mo(P.sat[rt_aval][rt_amont](tau)) * Transo_aval * phio_aval;
      
      R v =  P.Mo(X.S(o_amont)) * Transo_amont
	+ P.Mo(P.sat[rt_aval][rt_amont](tau)) * Transo_aval ;
      return u/v;
    } else {
      return phio_aval;
    }
  };

  std::function<R(R)> dPhio = [&](R tau) -> R { // & : capture all constants by ref
    if (P.Mo(X.S(o_amont)) > tolzero) {
      R u = P.Mo(X.S(o_amont)) * Transo_amont * phio_amont
	+ P.Mo(P.sat[rt_aval][rt_amont](tau)) * Transo_aval * phio_aval;
      
      R v =  P.Mo(X.S(o_amont)) * Transo_amont
	+ P.Mo(P.sat[rt_aval][rt_amont](tau)) * Transo_aval ;
      
      R du = P.dsat[rt_aval][rt_amont](tau) * P.dMo(P.sat[rt_aval][rt_amont](tau)) * Transo_aval * phio_aval ;
      R dv =  P.dsat[rt_aval][rt_amont](tau) * P.dMo(P.sat[rt_aval][rt_amont](tau)) * Transo_aval ;
      
      return (du*v - u*dv) / (v*v);
    } else {
      return 0.;
    }
  };
  
  // ------------------------------------------------------------------------------------- //
  // water
  R phiw_left = X.P(left)  + P.rhow * P.grav * Mh->Mailles[left].xT(1);
  R phiw_right = X.P(right) + P.rhow * P.grav * Mh->Mailles[right].xT(1);

  // maille amont/aval
  uint w_amont = (phiw_left >= phiw_right ? left : right);
  uint w_aval = (phiw_left >= phiw_right ? right : left);

  // rocktypes
  rt_amont = Mh->Mailles[w_amont].rt;
  rt_aval = Mh->Mailles[w_aval].rt;
  assert(rt_amont != rt_aval);
  
  // phi amont, aval
  R phiw_amont = X.P(w_amont)  + P.rhow * P.grav * Mh->Mailles[w_amont].xT(1);
  R phiw_aval  = X.P(w_aval)   + P.rhow * P.grav * Mh->Mailles[w_aval].xT(1);
  
  // Transmissibilities
  R Transw_amont = P.perm[rt_amont] / (F.xF - Mh->Mailles[w_amont].xT).norm();
  R Transw_aval = P.perm[rt_aval] / (F.xF - Mh->Mailles[w_aval].xT).norm();

  // Phi_w and its derivative
  std::function<R(R)> Phiw = [&](R tau) -> R { // = capture all constants by ref
    if (P.Mw(1-X.S(w_amont)) > tolzero) {
      R u = P.Mw(1-X.S(w_amont)) * Transw_amont * phiw_amont
	+ P.Mw(1 - P.sat[rt_aval][rt_amont](tau)) * Transw_aval * phiw_aval;
      
      R v =  P.Mw(1-X.S(w_amont)) * Transw_amont
	+ P.Mw(1 - P.sat[rt_aval][rt_amont](tau)) * Transw_aval ;
      
      return u/v;
    } else {
      return phiw_aval;
    }
  };
  
  std::function<R(R)> dPhiw = [&](R tau) -> R { // = : capture all constants by ref
    if (P.Mw(1-X.S(w_amont)) > tolzero) {
      R u = P.Mw(1-X.S(w_amont)) * Transw_amont * phiw_amont
	+ P.Mw(1-P.sat[rt_aval][rt_amont](tau)) * Transw_aval * phiw_aval;
      R du = - P.dsat[rt_aval][rt_amont](tau) * P.dMw(1-P.sat[rt_aval][rt_amont](tau)) * Transw_aval * phiw_aval ;
      
      R v =  P.Mw(1-X.S(w_amont)) * Transw_amont + P.Mw(1-P.sat[rt_aval][rt_amont](tau)) * Transw_aval ;
      R dv =   - P.dsat[rt_aval][rt_amont](tau) * P.dMw(1-P.sat[rt_aval][rt_amont](tau)) * Transw_aval ;
      return (du*v - u*dv) / (v*v);
    } else {
      return 0.;
    }
  };

  
  // -------------------------------------------------------------------- //
  //
  //     Solve
  //
  // -------------------------------------------------------------------- //
  const std::string logfile = "./output/newton1d_fails_n" + std::to_string(nt) + ".tmp";
  uint statusnwt = 1;
  uint statusdico = 1;
  
  // -------------------------------------------------------------------- //
  // Single-phase case : no gas 
  //
  if ( P.Mo(X.S(left)) + P.Mo(X.S(right)) < tolzero )
    {
      _tau = P.taumin; // defined in Physics
      _pw = Phiw(_tau) - P.rhow * P.grav * F.xF(1);
      _po = _pw + P.pc[rtl][rtr](_tau);

      F.istaumin = true;
      F.istaumax = false;
    }
  
  // -------------------------------------------------------------------- //  
  // Two-phase upwind degenerate - no oil in upwind
  //
  else if (P.Mo(X.S(o_amont)) < tolzero)
    {
      _tau = P.taumin;
      _pw = Phiw(_tau) - P.rhow * P.grav * F.xF(1);
      _po = _pw + P.pc[rtl][rtr](_tau);      
      R _phio = _po + P.rhoo * P.grav * F.xF(1);

      if ( _phio >= phio_aval ) {
	F.istaumin = true;
	F.istaumax = false;
      } else {

	F.istaumin = false;
	F.istaumax = false;
	
	// Residual functions
	R leftbound = P.taumin;
	R rightbound = P.taumax;
	
	std::function<R(R)> F_res = [&](R tau) -> R { return Phiw(tau) - P.rhow * P.grav * F.xF(1) + P.pc[rtl][rtr](tau) - phio_aval  + P.rhoo * P.grav * Mh->Mailles[o_aval].xT(1); };
	std::function<R(R)> dF_res = [&](R tau) -> R { return dPhiw(tau) + P.dpc[rtl][rtr](tau);  };
	statusdico = dichotomie(F_res, toldicho, leftbound, rightbound);	
	statusnwt = Newton1D(F_res, dF_res, leftbound, rightbound, tolnwt, nwt_max_iter, _tau, logfile);
      }
    }

  // -------------------------------------------------------------------- //  
  // Two-phase upwind degenerate - no oil in upwind
  //
  else if (P.Mo(X.S(o_aval)) < tolzero)
    {
      _tau = P.taumin;
      _pw = Phiw(_tau) - P.rhow * P.grav * F.xF(1);
      _po = _pw + P.pc[rtl][rtr](_tau);      
      R _phio = _po + P.rhoo * P.grav * F.xF(1);

      if ( _phio >= phio_amont ) {
	F.istaumin = true;
	F.istaumax = false;
      } else {

	F.istaumin = false;
	F.istaumax = false;
	
	// Residual functions
	R leftbound = P.taumin;
	R rightbound = P.taumax;
	
	std::function<R(R)> F_res = [&](R tau) -> R { return Phiw(tau) - P.rhow * P.grav * F.xF(1) + P.pc[rtl][rtr](tau) - phio_amont  + P.rhoo*P.grav*Mh->Mailles[o_amont].xT(1); };
	std::function<R(R)> dF_res = [&](R tau) -> R { return dPhiw(tau) + P.dpc[rtl][rtr](tau);  };

	statusdico = dichotomie(F_res, toldicho, leftbound, rightbound);	
	statusnwt = Newton1D(F_res, dF_res, leftbound, rightbound, tolnwt, nwt_max_iter, _tau, logfile);
      }
    }
  

  
  // -------------------------------------------------------------------- //
  // Two-phase non degenerate case
  else if ( P.Mo(X.S(o_amont)) > tolzero and P.Mo(X.S(o_aval)) > tolzero )
    {
      _tau = P.taumin;
      _pw = Phiw(_tau) - P.rhow * P.grav * F.xF(1);
      _po = _pw + P.pc[rtl][rtr](_tau);
      R _phio = _po + P.rhoo * P.grav * F.xF(1);

      if ( _phio >= phio_amont) {
	F.istaumin = true;
	F.istaumax = false;	
      } else {

	F.istaumin = false;
	F.istaumax = false;
	
	// Residual functions
	R leftbound = P.taumin;
	R rightbound = P.taumax;	
	
	std::function<R(R)> F_res = [&](R tau) -> R { return Phiw(tau) - P.rhow * P.grav * F.xF(1) + P.pc[rtl][rtr](tau) - Phio(tau)  + P.rhoo*P.grav*F.xF(1); };
	std::function<R(R)> dF_res = [&](R tau) -> R { return dPhiw(tau) + P.dpc[rtl][rtr](tau) - dPhio(tau); };

	statusdico = dichotomie(F_res, toldicho, leftbound, rightbound);	
	statusnwt = Newton1D(F_res, dF_res, leftbound, rightbound, tolnwt, nwt_max_iter, _tau, logfile);
      }

    }

  // -------------------------------------------------------------------- //
  else {
    std::cout << prefix  << "WARNING! This case should be treated : "
	      << "Mmo=" << P.Mo(X.S(o_aval))
      	      << " MMo=" << P.Mo(X.S(o_amont))
      	      << " Mmw=" << P.Mw(1-X.S(w_aval))
      	      << " MMw=" << P.Mw(1-X.S(w_amont))
	      << std::endl;
    exit(EXIT_FAILURE);
  }

  _pw = Phiw(_tau) - P.rhow * P.grav * F.xF(1);
  _po = _pw + P.pc[rtl][rtr](_tau);
  
  ResultVector << _tau, _po, _pw, P.pc[rtl][rtr](_tau);
  return ResultVector;
}
