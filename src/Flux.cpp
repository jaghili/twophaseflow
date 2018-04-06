#include "defs.h"
#include "Mesh.h"
#include "Physics.h"

R invtreshold (R x) {
  return (std::abs(x) < 1e-30 ? ((x>0.)-(x<0.))*1e-30 : 1./x); 
}

// Flux computation
void computeF(Mesh* Mh, Dat& X, Physics & P) {  // Compute Fluxes and Gradients
  // Inner fluxes
  uint iF = 0;
  for (std::vector<Face>::iterator F = Mh->FacesStd.begin(); F != Mh->FacesStd.end(); F++)
    { // loop over all inner faces

      uint left = F->iT_left;
      uint right = F->iT_right;
      Element* Tleft = &Mh->Mailles[F->iT_left];
      Element* Tright = &Mh->Mailles[F->iT_right];
      
      uint rtd = Tleft->rt; // lower rocktype
      uint rtu = Tright->rt; // upper rocktype

      R A = - ((Tleft->xT - Tright->xT).unaryExpr(&invtreshold).dot(F->nF));
      
      R s1 = P.perm[rtd];/// (F->xF - Tleft->xT).norm();
      R s2 = P.perm[rtu];// / (F->xF - Tright->xT).norm();

      R Trans = 2./(1./s1 + 1./s2) * A;

      F->flux_oil = Trans * ( P.rhoo/A * P.vg.dot(F->nF) + (X.P(left) + P.pc[rtd][rtd](X.S(left)) - X.P(right) - P.pc[rtu][rtu](X.S(right))) );
      F->flux_wat = Trans * ( P.rhow/A * P.vg.dot(F->nF) + X.P(left) - X.P(right) );
      iF++;
    }
  //std::cout << "# Loops in inner Fluxes:" << iF << std::endl;


  // Boundary fluxes
  int iFb=0;
  for(std::vector<BFace>::iterator F = Mh->BoundaryFaces.begin(); F != Mh->BoundaryFaces.end(); F++) {
    uint iT = F->iT_left;
    Element* T = &Mh->Mailles[iT];
    uint rt = T->rt; // get rocktype of bottom cell
    R A = -(T->xT - F->xF).unaryExpr(&invtreshold).dot(F->nF);
    R Trans = P.perm[rt] * A;// / (F->xF - T->xT).norm();

    F->flux_oil = Trans * (P.rhoo/A * P.vg.dot(F->nF) + (X.P(iT) + P.pc[rt][rt](X.S(iT)) - F->pw - P.pc[rt][rt](F->sat) )) ;    
    F->flux_wat = Trans * (P.rhow/A * P.vg.dot(F->nF) + X.P(iT) - F->pw) ;
    
    iFb++;
  }
  //std::cout << "# Loops on Bd Fluxes:" << iFb << std::endl;

  
  // Compute Fluxes on Interfaces
  uint iFi = 0;
  for (std::vector<IFace>::iterator F = Mh->Interfaces.begin(); F != Mh->Interfaces.end(); F++) {
    uint left = F->iT_left;
    uint right = F->iT_right;
    
    Element* Tleft = &Mh->Mailles[F->iT_left];
    Element* Tright = &Mh->Mailles[F->iT_right];
    
    uint rtd = Tleft->rt; // lower rocktype
    uint rtu = Tright->rt; // upper rocktype
    
    Vec rhog_minus_gradu = Vec::Zero(2);
    {// OIL
      // left
      R A = - ((Tleft->xT - F->xF).unaryExpr(&invtreshold).dot(F->nF));      
      R Trans_d = P.perm[rtd] * A; // / (F->xF - Tleft->xT).norm(); // left
      F->flux_oil_left = Trans_d * (P.rhoo/A * P.vg.dot(F->nF) + (X.P(left) + P.pc[rtd][rtd](X.S(left)) - X.Pi(iFi) - P.pc[rtd][rtu](X.Tau(iFi))) );

      // right
      A = (Tright->xT - F->xF).unaryExpr(&invtreshold).dot(F->nF);      
      R Trans_u = P.perm[rtu] * A; 
      F->flux_oil_right = Trans_u * (-P.rhoo/A * P.vg.dot(F->nF) + X.P(right) + P.pc[rtu][rtu](X.S(right)) - X.Pi(iFi) - P.pc[rtu][rtd](X.Tau(iFi)));
    }


    {// water
      //left
      R A = - ((Tleft->xT - F->xF).unaryExpr(&invtreshold).dot(F->nF));      
      R Trans_d = P.perm[rtd] * A; // / (F->xF - Tleft->xT).norm(); // left
      F->flux_wat_left = Trans_d * (P.rhow/A * P.vg.dot(F->nF) + (X.P(left) - X.Pi(iFi) ));
      //right
      A = ((Tright->xT - F->xF).unaryExpr(&invtreshold).dot(F->nF));      
      R Trans_u =  P.perm[rtu] * A; 
      F->flux_wat_right =  Trans_u * (-P.rhow/A * P.vg.dot(F->nF) + X.P(right) - X.Pi(iFi) );
    }
    iFi++;
  }
  //std::cout << "# Loops on interfaces Fluxes :" << iFi << std::endl;
} // End flux computation
