#include "BasicFunctions.h"
#include <fstream>
#include <iomanip>

// #define BF_DEBUG

// Dichotomie Solver
uint dichotomie(std::function<R(R)> & f, R tol, R & a, R & b) {
  char* prefix = "[Bisection] ";
  if (f(a)*f(b) > 0) {
    return 0;
  }
  else {
    // Find correct bounds
    int iter=0;

#ifdef BF_DEBUG
  std::cout.precision(16);
  std::cout << prefix
	    << std::left << std::setw(30) << std::scientific << "a"
	    << std::right << std::scientific << a
	    << std::endl;
  std::cout << prefix
	    << std::left << std::setw(30) << std::scientific << "b"
	    << std::right << std::scientific << b
	    << std::endl;
  std::cout << prefix
	    << std::left << std::setw(30) << std::scientific << "tol"
	    << std::right << std::scientific << tol
	    << std::endl;
#endif
  
    while (f(a)*f(b) < 0. and std::abs(f((a+b)/2)) > tol and iter < 1000) {
      iter++;
      if ( f((a+b)/2)*f(b) > 0 ) { b = (a+b)/2; }
      else { a = (a+b)/2; }
#ifdef BF_DEBUG
      std::cout << prefix
		<< std::setw(5) << std::scientific << iter
		<< std::setw(30) << std::scientific << a
		<< std::setw(30) << std::scientific << b
		<< std::setw(30) << std::scientific << f((a+b)/2.)
		<< std::endl;
#endif
    }
    return 1;
  }
}

// 1D Newton algorithm
// status = 1 (success) or 0 (fails)
uint Newton1D(std::function<R(R)> & f, // f
	      std::function<R(R)> & df, // dérivée
	      R xleft,
	      R xright,
	      R tol,
	      int maxiter,
	      R & sol,
	      std::string logfile) { // xinit
  
  uint status=1; // default = fail
  Vec shoots = Vec::Zero(maxiter);
  char* nwt = "[Newton1D] ";

#ifdef BF_DEBUG
  std::cout.precision(16);
  std::cout << nwt
	    << std::left << std::setw(30) << std::scientific << "a"
	    << std::right << std::scientific << xleft
	    << std::endl;
  std::cout << nwt
	    << std::left << std::setw(30) << std::scientific << "b"
	    << std::right << std::scientific << xright
	    << std::endl;
  std::cout << nwt
	    << std::left << std::setw(30) << std::scientific << "tol"
	    << std::right << std::scientific << tol
	    << std::endl;
 #endif
  
  R x = (xleft+xright)/2.;
  R err = std::abs(f(x)); 
  uint iter=0; // compteur
  R w = 1.;

#ifdef BF_DEBUG
  std::cout << nwt
	    << std::setw(5) << std::scientific << "iter"
	    << std::setw(30) << std::scientific << "x"
	    << std::setw(30) << std::scientific << "err"
	    << std::endl;
#endif
  while (err > tol and iter < maxiter) {
#ifdef BF_DEBUG
    std::cout << nwt
	      << std::setw(5) << std::scientific << iter
	      << std::setw(30) << std::scientific << x
      	      << std::setw(30) << std::scientific << err
	      << std::endl;
#endif
    R xold = x;
    x = xold - w * f(xold)/df(xold); // Newton update
    if (x<=xleft) { x = xleft+1e-10; }
    if (x>=xright) { x = xright-1e-10; }
    
    err = std::abs(x-xold);
    // nwt_file_ifacexver << nwt_offset_iface + iter << "\t" << err << std::endl;

    // update
    shoots(iter) = x;
    iter++;
    w = 0.95*w;
  }
#ifdef BF_DEBUG
  std::cout << nwt
	    << std::setw(5) << std::scientific << iter
	    << std::setw(30) << std::scientific << x
	    << std::setw(30) << std::scientific << err
    	    << std::endl;
#endif
  //nwt_file_ifacexver << "\n\n";
  //nwt_offset_iface += 3;
  
  std::cout.precision(5);
  // Newton has failed : writing in logfile
  if (iter == maxiter or err > tol) { 
    // writing log
    
    std::ofstream write_logfile(logfile, std::ofstream::app); // append file
    for (R _x=xleft; _x < xright; _x += (xright-xleft)/100.) {
      write_logfile << std::setprecision(15) << std::scientific
		    << _x << "\t" << f(_x) << "\t" << df(_x)
		    << std::endl;
    }
    write_logfile << "\n\n";
    
    for (int i = 0; i < maxiter; i++) {
      write_logfile << shoots(i) << "\t" << f(shoots(i)) << "\n"; // shoots
    }
    write_logfile << "\n\n\n";
    write_logfile.close();
    
    // error
    status = 0;
  }
  else {
    status = 1;
  }
  sol = x;
  return status;
}
