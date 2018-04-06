#include "Jacobian.h"
#include "defs.h"

extern int Ndof;
//extern int N;
//extern int Ni;
extern R dt;
extern R invtreshold(R x);
extern R vvol;

// Computation of the (Jacobian in sparse format
SpMat computeJacobian(Mesh* Mh, Dat & X, Dat & Xold, Physics & P) {
  //Mat dRes = Mat::Zero(2*Ndof,2*Ndof);
  int N = Mh->N;
  uint Ndof = X.Ndof;
  SpMat dRes(2*Ndof,2*Ndof);
  std::vector<Triplet> coeffs;
  coeffs.reserve(6 * 2 * Ndof); // estimation of nnz entries
  int iFi = 0;
  
  // ------------------------------------------------------------------------ //
  // Accumulation terms in cells
  for (std::vector<Element>::iterator M = Mh->Mailles.begin(); M != Mh->Mailles.end(); M++)
    {
      uint i = M->idx;
      uint oil = 2*i;
      uint wat = 2*i+1;
      
      uint dp = 2*i;
      uint ds = 2*i+1;
      
      // oil
      coeffs.push_back(Triplet(oil, dp, 0.));
      coeffs.push_back(Triplet(oil, ds, P.rhoo * P.phi * M->vol /dt));

      // water
      coeffs.push_back(Triplet(wat, dp, 0.));
      coeffs.push_back(Triplet(wat, ds, - P.rhow * P.phi * M->vol /dt));
    } // OK

  // Accumulation terms in interfaces (fictive) cells
  iFi = 0;
  for(std::vector<IFace>::iterator F = Mh->Interfaces.begin(); F != Mh->Interfaces.end(); F++)
    { // in all interfaces
      uint rtd = Mh->Mailles[F->iT_left].rt; // rt of lower cell
      uint rtu = Mh->Mailles[F->iT_right].rt; // rt of upper cell

      uint oil = 2*(N+iFi);
      uint wat = 2*(N+iFi)+1;
      uint dtau = 2*(N+iFi)+1;

      coeffs.push_back(Triplet(oil, dtau, P.rhoo * P.phi * vvol /dt * P.dsat[rtd][rtu](X.Tau(iFi)) + P.rhoo * P.phi * vvol /dt * P.dsat[rtu][rtd](X.Tau(iFi))));
      
      coeffs.push_back(Triplet(wat, dtau, - P.rhow * P.phi * vvol /dt * P.dsat[rtd][rtu](X.Tau(iFi)) -P.rhow * P.phi * vvol /dt * P.dsat[rtu][rtd](X.Tau(iFi))));
      iFi++;
    }
  
  
  // ------------------------------------------------------------------------ //
  // Adding fluxes : inner cells
  for (std::vector<Face>::iterator F = Mh->FacesStd.begin(); F != Mh->FacesStd.end(); F++)
    { // toutes les faces internes
      
      // indices
      uint left = F->iT_left;
      uint right = F->iT_right;

      Element* Tleft = &Mh->Mailles[left];
      Element* Tright = &Mh->Mailles[right];
      
      // rocktypes
      uint rtd = Tleft->rt; // lower rocktype
      uint rtu = Tright->rt; // upper rocktype

      // Same rocktypes
      if (rtd == rtu) { 
	R _s1 = P.perm[rtd];
	R _s2 = P.perm[rtu];

      	R A = - ((Tleft->xT - Tright->xT).unaryExpr(&invtreshold).dot(F->nF));            
	R Trans =  2./(1./_s1 + 1./_s2) * A ;

	uint oil_g = 2*left;
	uint oil_d = 2*right;
	uint wat_g = 2*left+1;
	uint wat_d = 2*right+1;

	uint dpdg = 2*left;
	uint dpdd = 2*right;
	uint dsdg = 2*left+1;
	uint dsdd = 2*right+1;

	uint jup1 = (F->flux_oil>0 ? left : right); // oil
	uint jup2 = (F->flux_wat>0 ? left : right); // wat

	uint dsdup1 = 2*jup1+1; // @s(jup1)
	uint dsdup2 = 2*jup2+1; // @s(jup2)
      
	// Derivative wrt p(i-1) et p(i)
	// Oil fluxes      
	R ss = Trans * P.rhoo * P.Mo(X.S(jup1));
	coeffs.push_back(Triplet(oil_g, dpdg, ss));
	coeffs.push_back(Triplet(oil_g, dpdd, -ss));
	coeffs.push_back(Triplet(oil_d, dpdg, -ss));
	coeffs.push_back(Triplet(oil_d, dpdd, ss));
      
	// water fluxes
	ss = Trans * P.rhow * P.Mw(1-X.S(jup2));
	coeffs.push_back(Triplet(wat_g, dpdg, ss));
	coeffs.push_back(Triplet(wat_g, dpdd, -ss));
	coeffs.push_back(Triplet(wat_d, dpdg, -ss));
	coeffs.push_back(Triplet(wat_d, dpdd, ss));

	// Derivatives wrt S
	// oil eqns
	R s1 = Trans * P.rhoo * P.Mo(X.S(jup1)) * P.dpc[rtd][rtd](X.S(left));
	R s2 = Trans * P.rhoo * P.Mo(X.S(jup1)) * P.dpc[rtu][rtu](X.S(right));
	coeffs.push_back(Triplet(oil_g, dsdg, s1));
	coeffs.push_back(Triplet(oil_g, dsdd, -s2));
	coeffs.push_back(Triplet(oil_d, dsdg, -s1));
	coeffs.push_back(Triplet(oil_d, dsdd, s2));
      
	// derivative wrt to S upwind
	//oil
	ss = F->flux_oil * P.rhoo * P.dMo(X.S(jup1));
	coeffs.push_back(Triplet(oil_g, dsdup1, ss));
	coeffs.push_back(Triplet(oil_d, dsdup1, -ss));
      
	//water
	ss = F->flux_wat * P.rhow * P.dMw(1-X.S(jup2));
	coeffs.push_back(Triplet(wat_g, dsdup2, -ss));
	coeffs.push_back(Triplet(wat_d, dsdup2, ss));      
      }
    }
  
  // Toutes les interfaces
  iFi = 0;
  for (std::vector<IFace>::iterator F = Mh->Interfaces.begin(); F != Mh->Interfaces.end(); F++) {
    // ----------------------------------------------------------------------- //    

    uint rtd = Mh->Mailles[F->iT_left].rt;
    uint rtu = Mh->Mailles[F->iT_right].rt;
      
    {// interface gauche entre maille et mini-maille 
      uint left = F->iT_left;
      uint right = N + iFi;
      
      uint oil_g = 2*left;
      uint oil_d = 2*right;
      uint wat_g = 2*left+1;
      uint wat_d = 2*right+1;

      uint dpdg = 2*left;
      uint dsdg = 2*left+1;
      uint dpi = 2*right;
      uint dtau = 2*right+1;

      Element* Tleft = &Mh->Mailles[left];
      
      R A = - ((Tleft->xT - F->xF).unaryExpr(&invtreshold).dot(F->nF));      
      R Trans = P.perm[rtd] * A;
      
      R mob, dmob; // mobility and its derivative
      uint jup; // upwind element

      // oil
      if (F->flux_oil_left >= 0.) {
	mob = P.rhoo * P.Mo(X.S(left));
	dmob = P.rhoo * P.dMo(X.S(left)); // dmob/dsL
	jup = left;
      } else {
	mob = P.rhoo * P.Mo(P.sat[rtd][rtu](X.Tau(iFi)));
	dmob = P.rhoo * P.dMo(P.sat[rtd][rtu](X.Tau(iFi))) * P.dsat[rtd][rtu](X.Tau(iFi));
	jup = right;
      }

      // dérivées par rapport à p(left) et pi
      R ss = Trans * mob;
      coeffs.push_back(Triplet(oil_g, dpdg, ss));
      coeffs.push_back(Triplet(oil_g, dpi, -ss));
      coeffs.push_back(Triplet(oil_d, dpdg, -ss));
      coeffs.push_back(Triplet(oil_d, dpi, ss));

	
      // dérivées de la mobilité mob1 par rapport à S ou Tau 
      ss = dmob * F->flux_oil_left;
      coeffs.push_back(Triplet(oil_g, 2*jup+1, ss));
      coeffs.push_back(Triplet(oil_d, 2*jup+1, -ss));

	
      // dérivées par rapport à Pc(S(left)) par rapport à S  
      ss = Trans * mob * P.dpc[rtd][rtd](X.S(left));
      coeffs.push_back(Triplet(oil_g, dsdg, ss));
      coeffs.push_back(Triplet(oil_d, dsdg, -ss));

      // dérivées part rapport à pc[rtd][rtu](tau(iFi))
      ss = - Trans * mob * P.dpc[rtd][rtu](X.Tau(iFi));
      coeffs.push_back(Triplet(oil_g, dtau, ss));
      coeffs.push_back(Triplet(oil_d, dtau, -ss));

      // water
      if (F->flux_wat_left >= 0.) {
	mob = P.rhow * P.Mw(1.-X.S(left));
	dmob = - P.rhow * P.dMw(1.-X.S(left));
	jup = left;
      } else {
	mob = P.rhow * P.Mw(1.-P.sat[rtd][rtu](X.Tau(iFi)));
	dmob = - P.rhow * P.dMw(1.-P.sat[rtd][rtu](X.Tau(iFi))) * P.dsat[rtd][rtu](X.Tau(iFi));
	jup = right;
      }

      // dérivées par rapport à p(left) et pi
      ss = Trans * mob;
      coeffs.push_back(Triplet(wat_g, dpdg, ss));
      coeffs.push_back(Triplet(wat_g, dpi, -ss));
      coeffs.push_back(Triplet(wat_d, dpdg, -ss));
      coeffs.push_back(Triplet(wat_d, dpi, ss));
	
      // dérivées par rapport à la mobilité
      ss = dmob * F->flux_wat_left;
      coeffs.push_back(Triplet(wat_g, 2*jup+1, ss));
      coeffs.push_back(Triplet(wat_d, 2*jup+1, -ss));
    }
    // ---------------------------------------------------------------------------------- //
    {       // interface entre minimaille et maille droite, oil
      uint left = N + iFi;
      uint right = F->iT_right;
      
      uint oil_g = 2*left;
      uint oil_d = 2*right;
      uint wat_g = 2*left+1;
      uint wat_d = 2*right+1;

      uint dpdd = 2*right;
      uint dsdd = 2*right+1;
      uint dpi  = 2*left;
      uint dtau = 2*left+1;

      Element* Tright = &Mh->Mailles[right];
      
      R A = ((Tright->xT - F->xF).unaryExpr(&invtreshold).dot(F->nF));            
      R Trans = P.perm[rtu] * A;

      // oil
      R mob, dmob; // mobility and its derivative
      uint jup; // upwind element

      if (F->flux_oil_right >= 0.) {
	mob = P.rhoo * P.Mo(X.S(right));
	dmob = P.rhoo * P.dMo(X.S(right));
	jup = right; 
      } else {
	mob  = P.rhoo * P.Mo(P.sat[rtu][rtd](X.Tau(iFi)));
	dmob = P.rhoo * P.dMo(P.sat[rtu][rtd](X.Tau(iFi))) * P.dsat[rtu][rtd](X.Tau(iFi)); // dmob/dtau
	jup = left;
      }

      // dérivées par rapport à p(right) et pi
      R ss = Trans * mob;
      coeffs.push_back(Triplet(oil_d, dpdd, ss));
      coeffs.push_back(Triplet(oil_d, dpi, -ss));
      coeffs.push_back(Triplet(oil_g, dpdd, -ss));
      coeffs.push_back(Triplet(oil_g, dpi, ss));
      
      // dérivées de la mobilité par rapport à S
      ss = dmob * F->flux_oil_right;
      coeffs.push_back(Triplet(oil_d, 2*jup+1, ss));
      coeffs.push_back(Triplet(oil_g, 2*jup+1, -ss));
	
      // dérivées par rapport à Pc(S(right))
      ss = Trans * mob * P.dpc[rtu][rtu](X.S(right));
      coeffs.push_back(Triplet(oil_d, dsdd, ss));
      coeffs.push_back(Triplet(oil_g, dsdd, -ss));
      
	
      // dérivées par rapport à Pc[rtu][rtd]
      ss = - Trans *  mob * P.dpc[rtu][rtd](X.Tau(iFi));
      coeffs.push_back(Triplet(oil_d, dtau, ss));
      coeffs.push_back(Triplet(oil_g, dtau, -ss));

      
      // water
      if (F->flux_wat_right >= 0.) {
	mob = P.rhow * P.Mw(1.-X.S(right));
	dmob = - P.rhow * P.dMw(1.-X.S(right));
	jup = right;
      } else {
	mob = P.rhow * P.Mw(1.-P.sat[rtu][rtd](X.Tau(iFi)));
	dmob = - P.rhow * P.dMw(1.-P.sat[rtu][rtd](X.Tau(iFi))) * P.dsat[rtu][rtd](X.Tau(iFi));
	jup = left;
      }
      
      // dérivées par rapport à pj(right) et pi      
      ss = Trans * mob;
      coeffs.push_back(Triplet(wat_d,dpdd, ss));
      coeffs.push_back(Triplet(wat_d,dpi, -ss));
      coeffs.push_back(Triplet(wat_g,dpdd,-ss));
      coeffs.push_back(Triplet(wat_g,dpi,  ss));
      
      // dérivées par rapport à la mobilité
      ss = dmob * F->flux_wat_right;
      coeffs.push_back(Triplet(wat_d,2*jup+1, ss));
      coeffs.push_back(Triplet(wat_g,2*jup+1,-ss));
    }
      
    // end
    iFi++;
    } // end loop interfaces 
  
  // ---------------------------
  // boundary terms
  for (std::vector<BFace>::iterator F = Mh->BoundaryFaces.begin(); F != Mh->BoundaryFaces.end(); F++) {
    if (F->DIRICHLET_BC) {
      Element* T = &Mh->Mailles[F->iT_left];
      uint ib = T->idx;
      uint rt = T->rt; // get rocktype of bottom cell

      R A = -(T->xT - F->xF).unaryExpr(&invtreshold).dot(F->nF);
      R Trans = P.perm[rt] * A;
      
      uint oil = 2*ib;
      uint wat = 2*ib+1;
      uint dp = 2*ib;
      uint ds = 2*ib+1;
      
      // oil
      if (F->flux_oil < 0.) {
	coeffs.push_back( Triplet(oil, dp,  Trans * P.rhoo * P.Mo(F->sat)) );
	coeffs.push_back( Triplet(oil, ds,  Trans * P.rhoo * P.Mo(F->sat) * P.dpc[rt][rt](X.S(ib))) );
      } else {
	coeffs.push_back( Triplet(oil, dp,  Trans * P.rhoo * P.Mo(X.S(ib)) ));
	coeffs.push_back( Triplet(oil, ds,  P.rhoo * P.dMo(X.S(ib)) * F->flux_oil ));
	coeffs.push_back( Triplet(oil, ds,  Trans * P.rhoo * P.Mo(X.S(ib)) * P.dpc[rt][rt](X.S(ib)) ));
      }
      // water
      if(F->flux_wat < 0) {
	coeffs.push_back( Triplet(wat, dp, Trans * P.rhow * P.Mw(1.- F->sat) ));
      } else {
	coeffs.push_back( Triplet(wat, dp,   Trans * P.rhow * P.Mw(1.-X.S(ib))));
	coeffs.push_back( Triplet(wat, ds,   - P.rhow * P.dMw(1.- X.S(ib)) * F->flux_wat ));
      }
    }    
  }
  
  //assemble dRes
  dRes.setFromTriplets(coeffs.begin(), coeffs.end());

  // display
  /*
  std::cout << "Jacobian=\n";
  for (int k=0; k<dRes.outerSize(); ++k) {
    for (SpMat::InnerIterator it(dRes,k); it; ++it)
      {
	std::cout << it.row() << "\t" << it.col() << "\t"<< it.value() << "\n";
      }
    std::cout << "\n";
  }
  */  
  return dRes;
}


// Computation of the (finite difference) Jacobian in dense format
Mat computeApproxJacobian (Mesh* Mh, Dat & X, Dat & Xold, Physics & P) {
  R h = 1e-9; // bug quand 1e-11
  uint Ndof = Mh->N + Mh->Ni;
  Mat Id = Mat::Identity(2*Ndof,2*Ndof);
  Mat dRes = Mat::Zero(2*Ndof, 2*Ndof);

  // Dat X1,X2;
  // for (int j = 0; j < dRes.cols() ; j++) {
  //   Vec vX2 = X.get() + h*Id.col(j);
  //   Vec vX1 = X.get(); 
  //   X2.set(vX2);
  //   X1.set(vX1);
  //   F2.computeF(Mh, X2, P);
  //   F1.computeF(Mh, X1, P);
  //   Mat _dRes  = ( computeResidual(Mh, F2, X2, Xold, P) - computeResidual(Mh, F1, X1, Xold, P) ) / h;
  //   for (uint i = 0; i < Ndof; i++) {
  //     dRes(2*i,j) = _dRes(i,0);
  //     dRes(2*i+1,j) = _dRes(i,1);
  //   }
  // }
  return dRes;
}
