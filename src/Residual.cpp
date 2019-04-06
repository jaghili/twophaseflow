#include "Residual.h"

extern R dt; // dt is defined outside of this file
extern R vvol;

// Compute Residual *onces fluxes are computed*
Mat computeResidual(Mesh* Mh, Dat & X, Dat & Xold, Physics & P) {
  assert(X.Ndof == Xold.Ndof);
  Mat Res = Mat::Zero(X.Ndof,2);

  // --------------------------------------------------------------------- //
  // Accumulation terms in all cells
  for (std::vector<Element>::iterator M = Mh->Mailles.begin(); M != Mh->Mailles.end(); M++) // loop on cells
    { 
      uint i = M->idx;
      Res(i,0) = P.rhoo * P.phi[M->rt] * (M->vol)/dt  * (X.S(i) - Xold.S(i)); // oil saturation
      Res(i,1) = - P.rhow * P.phi[M->rt] * M->vol/dt * (X.S(i) - Xold.S(i)); // water saturation

      /*std::cout << "Acc Elem no." << i
	<< "\t oil " << Res(i,0)
	<< "\t wat " << Res(i,1) << std::endl; */
    }

  // in all interfaces
  int iFi=0;
  for (std::vector<IFace>::iterator F = Mh->Interfaces.begin(); F != Mh->Interfaces.end(); F++) // loop on interfaces
    { 
      uint rtd = Mh->Mailles[F->iT_left].rt; // rt of lower cell
      uint rtu = Mh->Mailles[F->iT_right].rt; // rt of upper cell

      Res(X.N + iFi,0) =
	P.rhoo*P.phi[rtd]*vvol/ dt * (P.sat[rtd][rtu](X.Tau(iFi))-P.sat[rtd][rtu](Xold.Tau(iFi)))
	+P.rhoo*P.phi[rtu]*vvol/ dt * (P.sat[rtu][rtd](X.Tau(iFi))- P.sat[rtu][rtd](Xold.Tau(iFi)));

      Res(X.N + iFi,1) =
	-P.rhow *P.phi[rtd]*vvol/dt*(P.sat[rtd][rtu](X.Tau(iFi))-P.sat[rtd][rtu](Xold.Tau(iFi)))
	-P.rhow *P.phi[rtu]*vvol/dt*(P.sat[rtu][rtd](X.Tau(iFi))-P.sat[rtu][rtd](Xold.Tau(iFi)));
      
      /*std::cout << "Acc IFace no." << iFi
		<< "\t oil " << Res(X.N+iFi,0)
		<< "\t wat " << Res(X.N+iFi,1) << std::endl;*/
      
      iFi++;
    }
  
  // --------------------------------------------------------------------- //    
  // Décentrages on Standard Faces
  iFi = 0;
  for (std::vector<Face>::iterator F = Mh->FacesStd.begin(); F != Mh->FacesStd.end(); F++) // loop on std. faces
    { // toutes les faces internes
      
      uint left = F->iT_left;
      uint right = F->iT_right;
      
      uint rtd = Mh->Mailles[left].rt; // lower rocktype
      uint rtu = Mh->Mailles[right].rt; // upper rocktype
      
      if (rtd == rtu) {
	// Adding oil fluxes
	uint jup1 = (F->flux_oil >=0. ? left : right);
	R s1 = F->flux_oil * P.rhoo * P.Mo(X.S(jup1)); // oil flux
	Res(left,0) += s1;
	Res(right,0) -= s1;
      
	// Adding water fluxes
	uint jup2 = (F->flux_wat >=0. ? left : right);
	R s2 = F->flux_wat * P.rhow * P.Mw(1-X.S(jup2));
	Res(left,1) += s2;
	Res(right,1) -= s2;

	/*std::cout << "Dec StdFace no." << iFi
		<< "\t oil " << Res(left,0)
		<< "\t wat " << Res(right,1) << std::endl;*/
	
	iFi++;
      } 
    }
  
  // Décentrages from Interfaces
  iFi =0;
  for (std::vector<IFace>::iterator F = Mh->Interfaces.begin(); F!= Mh->Interfaces.end(); F++)  // loop on interfaces
    {
      uint left = F->iT_left;
      uint right = F->iT_right;
    
      uint rtd = Mh->Mailles[left].rt; // lower rocktype
      uint rtu = Mh->Mailles[right].rt; // upper rocktype

      // gauche
      R Fo_i_g = max(F->flux_oil_left, 0.) * P.rhoo * P.Mo(X.S(left)) + min(F->flux_oil_left,0.) * P.rhoo * P.Mo(P.sat[rtd][rtu](X.Tau(iFi)));
      Res(left,0) += Fo_i_g;
      if (!F->istaumin) { Res(X.N + iFi,0) -= Fo_i_g; } else { Res(X.N+iFi,0) = X.Tau(iFi) - P.taumin; }


      R Fw_i_g = max(F->flux_wat_left,0.) * P.rhow * P.Mw(1-X.S(left)) + min(F->flux_wat_left,0.) * P.rhow * P.Mw(1-P.sat[rtd][rtu](X.Tau(iFi))); 
      Res(left, 1) += Fw_i_g;
      Res(X.N + iFi, 1) -= Fw_i_g;
      
      // droite
      R Fo_i_d = min(F->flux_oil_right, 0.) * P.rhoo * P.Mo(P.sat[rtu][rtd](X.Tau(iFi))) + max(F->flux_oil_right,0.) * P.rhoo * P.Mo(X.S(right));
      Res(right,0) += Fo_i_d;
      if (!F->istaumin) { Res(X.N + iFi,0) -= Fo_i_d; } else { Res(X.N+iFi,0) = X.Tau(iFi) - P.taumin; }

      R Fw_i_d = min(F->flux_wat_right,0.) * P.rhow * P.Mw(1-P.sat[rtu][rtd](X.Tau(iFi))) + max(F->flux_wat_right,0.) * P.rhow * P.Mw(1-X.S(right)); // 
      Res(right,1) += Fw_i_d;
      Res(X.N + iFi,1) -= Fw_i_d;

      
      iFi++;
    }

  // Déc, on boundary faces
  for(std::vector<BFace>::iterator F = Mh->BoundaryFaces.begin(); F != Mh->BoundaryFaces.end(); F++)
    {
      uint iT = F->iT_left;
      if (F->DIRICHLET_BC) {
	/*
	std::cout << "*** CHECK ***" << std::endl;
	std::cout << "flux oil " << F->flux_oil << std::endl;
	std::cout << "rhoo " << P.rhoo << std::endl;
	std::cout << "Mo(satiF)" << F->sat << "\t" << P.Mo(F->sat) << std::endl;
	std::cout << "Mo(satC) " << X.S(iT) <<  std::endl;
	*/
	
	R s1 = min(F->flux_oil,0.) * P.rhoo * P.Mo(F->sat) + max(F->flux_oil,0.) * P.rhoo * P.Mo(X.S(iT));
	R s2 = min(F->flux_wat,0.) * P.rhow * P.Mw(1-F->sat) + max(F->flux_wat,0.) * P.rhow * P.Mw(1-X.S(iT));
	Res(iT,0) += s1; 
	Res(iT,1) += s2;

	/* std::cout << "Bound Face no." << iT
		<< "\t oil " << Res(iT,0)
		<< "\t wat " << Res(iT,1) << std::endl;	*/
	
      }
      /*
      // to do
      else if (F->HAS_NEUMANN_BC) {
	Res(iT,0) += min(F->neumann_flux_oil,0.) * P.rhoo * P.Mo(F->sat) + max(F->neumann_flux_oil,0.) * P.rhoo * P.Mo(X.S(iT)); //bottom oil
	Res(iT,1) += min(F->neumann_flux_wat,0.) * P.rhow * P.Mw(1-F->sat) + max(F->neumann_flux_wat,0.) * P.rhow * P.Mw(1-X.S(iT));
      }
      */
    }
  return Res;
} 
