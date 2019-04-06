#include "Mesh.h"

// compute Interfaces
void Mesh::computeInterfaces() {
  assert(ELEMENTS_ARE_SORTED);
  Ni=0; // interfaces
  
  // fill interfaces
  for (std::vector<Face>::iterator F = FacesStd.begin(); F != FacesStd.end(); F++) {
    if (Mailles[F->iT_left].rt != Mailles[F->iT_right].rt) { 

      IFace iF;
      iF.idx = F->idx;
      iF.iT_left = F->iT_left;
      iF.iT_right = F->iT_right;
      
      iF.xF = Vec::Zero(d);
      iF.xF(0)= F->xF(0);
      iF.xF(1)= F->xF(1);

      iF.xa = F->xa;
      iF.xb = F->xb;
      
      iF.vol = F->vol;
      iF.nF = Vec::Zero(d);
      iF.nF(0) = F->nF(0);
      iF.nF(1) = F->nF(1);
      Interfaces.push_back(iF); // ajouter F aux interfaces
      Ni++; // incrémenter
    }
  }
  //std::cout << "Number of interfaces: " << Ni << std::endl;
}

//extern int dim;
void Mesh::setRocktype(uint _idx, uint rocktype) {
  assert(_idx < N);
  assert(rocktype >= 0 && rocktype < Nrt);

  std::vector<Element>::iterator M = Mailles.begin();
  while (M->idx != _idx) { M++; }
  M->rt = rocktype;  
}


// display in terminal
void Mesh::display() {

  // Elements
  std::cout << LONGLINE << std::endl;
  for (auto & T : Mailles) {
    std::cout << "Element \tidx=" << T.idx << "\trt=" << T.rt << "\n";
  }
  for (auto & F : FacesStd) {
    std::cout << "Inner Face \tidx=" << F.idx
	      << " \tTleft=" << F.iT_left << "\tTright=" << F.iT_right
	      << "\tflux o=" << F.flux_oil << " w=" << F.flux_wat
	      << "\tnormal=" << F.nF(0) << " " << F.nF(1)
	      << "\n";
  }
  for (auto & F : BoundaryFaces) {
    std::cout << "Bd Face \tidx=" << F.idx << "\tT=" << F.iT_left;
    if (F.DIRICHLET_BC) { std::cout << "\tDIRICHLET BC  s=" << F.sat << " p=" << F.pw; }
    else if (F.NEUMANN_BC) { std::cout << "\tNEUMANN BC"; }
    std::cout << "\tflux o=" << F.flux_oil << " w=" << F.flux_wat
      	      << "\tnormal=" << F.nF(0) << " " << F.nF(1)
	      <<"\n";
  }
  for (auto & F : Interfaces) {
    std::cout << "Interface \tid=" << F.idx << " \tTleft=" << F.iT_left << "\tTright=" << F.iT_right
	      << "\tleft flux o=" << F.flux_oil_left << " w=" << F.flux_wat_left
      	      << " right flux o=" << F.flux_oil_right << " w=" << F.flux_wat_right
      	      << "\tnormal=" << F.nF(0) << " " << F.nF(1)
	      << "\n";
  }
  std::cout << LONGLINE << std::endl;
}
