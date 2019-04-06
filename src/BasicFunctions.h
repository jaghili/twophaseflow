#include "defs.h"
#include <functional>

// stop dichotomy when f((a+b)/2) < tol
uint dichotomie(std::function<R(R)> & f, R tol, R & a, R & b);

// 1D Newton algorithm
// xleft and xright are replaced by the final solution x
uint Newton1D(std::function<R(R)> & f,
	      std::function<R(R)> & df,
	      R xleft,
	      R xright,
	      R tol,
	      int maxiter,
	      R & sol,
	      std::string logfile);
