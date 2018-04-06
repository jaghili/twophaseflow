#include <iomanip>
#include "Data.h"
#include "Mesh.h"

Dat::Dat(Mesh* Mh) {
  N = Mh->N;
  Ni = Mh->Ni;

  Ndof = (N+Ni);
  
  S = Vec::Zero(N);
  P = Vec::Zero(N);
  Pi = Vec::Zero(Ni);
  Tau = Vec::Zero(Ni);
};

Vec Dat::get() {
  Vec U = Vec::Zero(2*N + 2*Ni);
  for (uint i = 0 ; i < N; i++) {
    U(2*i) = P(i);
    U(2*i+1) = S(i);      
  }
  for (uint j = 0; j < Ni; j++) {
    U(2 * N + 2*j ) = Pi(j);
    U(2 * N + 2*j + 1) = Tau(j);
      
  }
  return U;
};

void Dat::set(Vec & U) {
  for (uint i = 0 ; i < N; i++) {
    P(i) = U(2*i);
    S(i) = U(2*i+1);      
  }
  for (uint j = 0; j < Ni; j++) {
    Pi(j) = U(2*N + 2*j);
    Tau(j) = U(2*N + 2*j+1);
  }
};

void Dat::display() {
  for (uint i=0; i < N ; i++) {
    std::cout << std::setprecision(3)
	      << std::scientific
	      << S(i) << "\t"
	      << P(i) << "\t"
	      << (i<Ni?Pi(i):0.) << "\t"
	      << (i<Ni?Tau(i):0.) << std::endl;
  }
}
