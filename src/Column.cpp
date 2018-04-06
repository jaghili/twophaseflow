#include "Face.h"
#include "BFace.h"
#include "IFace.h"
#include "Element.h"
#include "Column.h"
#include <algorithm>

Column::Column(int n, double w, Physics & P) {
  assert(n>1);
  assert(P.a != P.b);
  assert(w > 0.);
  
  N = n; // elements
  Nf = n-1; // standard faces
  Nfb = n+n+2; // boundary elements
  
  // rocktypes info
  Ni = 0; // no. interfaces
  Nrt = 2;
  
  // alloc
  FacesStd.reserve(Nf);
  Mailles.reserve(N); // cells
  Interfaces.reserve(Nf); // Maximum number of interfaces
  BoundaryFaces.reserve(Nfb); // Boundary faces
  
  // Standard Faces
  R h = (P.b-P.a)/N;

  d = P.dim;
  
  // Faces standards
  for (uint i = 0; i < Nf; i++) {
    Face F;
    // face index
    F.idx = i;
    // neighboors
    F.iT_left = i;
    F.iT_right = i+1;
    // Face center
    F.xF = Vec::Zero(d);
    F.xF(0)= 0.5;
    F.xF(1)= P.a + h + i*h;
    F.xa = Vec::Zero(d);
    F.xb = Vec::Zero(d);
    F.xa(0)= 0.;
    F.xa(1)= P.a + h + i*h;
    F.xb(0)= w;
    F.xb(1)= P.a + h + i*h;
    // face volume
    F.vol = w;
    // normal
    F.nF = Vec::Zero(d);
    F.nF(0) = 0;
    F.nF(1) = 1;
    // add face
    FacesStd.push_back(F);
  }
  
  // fill mailles
  // inner elements
  for (uint i = 1; i < N-1; i++) {
    Element T;
    T.idx = i;
    T.xT = Vec::Zero(d);
    T.xT(0)=0.5;
    T.xT(1)=P.a + h/2. + i*h;
    T.rt = 0;
    T.vol = h*w;
    T.faces = iVec::Zero(4);
    T.faces(0) = i-1; // down
    T.faces(1) = i; // up
    T.faces(2) = n+1+i; // left
    T.faces(3) = 2*n+1+i; // right
    T.shape.setPointCount(4);
    T.shape.setPoint(0, sf::Vector2f(0., P.L - P.a - i*h ));
    T.shape.setPoint(1, sf::Vector2f(w, P.L - P.a - i*h ));
    T.shape.setPoint(2, sf::Vector2f(w, P.L - P.a - (i+1)*h ));
    T.shape.setPoint(3, sf::Vector2f(0., P.L - P.a - (i+1)*h ));
    
    // left boundary
    BFace F_l;
    F_l.idx = N+1+i;
    F_l.iT_left = i;
    F_l.iT_right = i;
    F_l.xF = Vec::Zero(d);
    F_l.xF(0)= 0.;
    F_l.xF(1)= P.a + h/2. + i*h;
    F_l.xa = Vec::Zero(d);
    F_l.xb = Vec::Zero(d);
    F_l.xa(0)= 0.;
    F_l.xa(1)= P.a + i*h;
    F_l.xb(0)= 0.;
    F_l.xb(1)= P.a + i*h + h;
    F_l.vol = h;
    F_l.nF = Vec::Zero(d);
    F_l.nF(0) = -1;
    F_l.nF(1) = 0;
    F_l.NEUMANN_BC = true;

    // right boundary
    BFace F_r;
    F_r.idx = 2*N+1+i;
    F_r.iT_left = i;
    F_r.iT_right = i;
    F_r.xF = Vec::Zero(d);
    F_r.xF(0)= 1;
    F_r.xF(1)= P.a + h/2. + i*h;
    F_r.xa = Vec::Zero(d);
    F_r.xb = Vec::Zero(d);
    F_r.xa(0)= w;
    F_r.xa(1)= P.a + i*h;
    F_r.xb(0)= w;
    F_r.xb(1)= P.a + i*h + h;
    F_r.vol = h;
    F_r.nF = Vec::Zero(d);
    F_r.nF(0) = 1.;
    F_r.nF(1) = 0.;
    F_r.NEUMANN_BC = true;
    
    Mailles.push_back(T);
    BoundaryFaces.push_back(F_l);
    BoundaryFaces.push_back(F_r);
  }

  // top/bottom elements
  {
    Element T_bottom;
    T_bottom.idx = 0;
    T_bottom.xT = Vec::Zero(d);
    T_bottom.xT(0)=0.5;
    T_bottom.xT(1)=P.a + h/2.;
    T_bottom.rt = 0;
    T_bottom.vol = h*w;
    T_bottom.faces = iVec::Zero(4);
    T_bottom.faces(0) = Nfb;
    T_bottom.faces(1) = 0;
    T_bottom.faces(2) = N+1; //left
    T_bottom.faces(3) = 2*N+1;
    T_bottom.shape.setPointCount(4);
    T_bottom.shape.setPoint(0, sf::Vector2f(0., P.L - P.a ));
    T_bottom.shape.setPoint(1, sf::Vector2f(w, P.L - P.a ));
    T_bottom.shape.setPoint(2, sf::Vector2f(w, P.L - P.a - h ));
    T_bottom.shape.setPoint(3, sf::Vector2f(0., P.L - P.a - h ));
    // bottom
    BFace F_bottom;
    F_bottom.idx = N-1;
    F_bottom.iT_left = 0;
    F_bottom.iT_right = 0;
    F_bottom.xF = Vec::Zero(d);
    F_bottom.xF(0) = 0.5;
    F_bottom.xF(1) = P.a;
    F_bottom.xa = Vec::Zero(d);
    F_bottom.xb = Vec::Zero(d);
    F_bottom.xa(0)= P.a;
    F_bottom.xa(1)= 0.;
    F_bottom.xb(0)= w;
    F_bottom.xb(1)= 0.;
    F_bottom.vol = 1.;
    F_bottom.nF = Vec::Zero(d);
    F_bottom.nF(0) = 0;
    F_bottom.nF(1) = -1.;
    F_bottom.pw = P.pw_bottom;
    F_bottom.sat = P.s_bottom;
    F_bottom.DIRICHLET_BC = true;    
    
    BFace F_bottom_l;
    F_bottom_l.idx = N+1;
    F_bottom_l.iT_left = 0;
    F_bottom_l.iT_right = 0;
    F_bottom_l.xF = Vec::Zero(d);
    F_bottom_l.xF(0)= 0.;
    F_bottom_l.xF(1)= P.a + h/2.;
    F_bottom_l.xa = Vec::Zero(d);
    F_bottom_l.xb = Vec::Zero(d);
    F_bottom_l.xa(0)= P.a;
    F_bottom_l.xa(1)= 0.;
    F_bottom_l.xb(0)= P.a;
    F_bottom_l.xb(1)= h;
    F_bottom_l.vol = h;
    F_bottom_l.nF = Vec::Zero(d);
    F_bottom_l.nF(0) = -1;
    F_bottom_l.nF(1) = 0;
    F_bottom.NEUMANN_BC = true;
    
    BFace F_bottom_r;
    F_bottom_r.idx = 2*N+1;
    F_bottom_r.iT_left = 0;
    F_bottom_r.iT_right = 0;
    F_bottom_r.xF = Vec::Zero(d);
    F_bottom_r.xF(0)= 1;
    F_bottom_r.xF(1)= P.a + h/2;
    F_bottom_r.xa = Vec::Zero(d);
    F_bottom_r.xb = Vec::Zero(d);
    F_bottom_r.xa(0)= w;
    F_bottom_r.xa(1)= 0.;
    F_bottom_r.xb(0)= w;
    F_bottom_r.xb(1)= P.a + h;
    F_bottom_r.vol = h;
    F_bottom_r.nF = Vec::Zero(d);
    F_bottom_r.nF(0) = 1.;
    F_bottom_r.nF(1) = 0;
    F_bottom.NEUMANN_BC = true;    
    
    // up
    Element T_up;
    T_up.idx = N-1;
    T_up.xT = Vec::Zero(d);
    T_up.xT(0)=0.5;
    T_up.xT(1)= P.b - h/2;
    T_up.rt = 0;
    T_up.vol = h*w;
    T_up.faces = iVec::Zero(4);
    T_up.faces(0) = N-2;
    T_up.faces(1) = Nfb + 1;
    T_up.faces(2) = 2*N;
    T_up.faces(3) = 3*N;
    T_up.shape.setPointCount(4);
    T_up.shape.setPoint(0, sf::Vector2f(0., P.L - P.b ));
    T_up.shape.setPoint(1, sf::Vector2f(w, P.L - P.b ));
    T_up.shape.setPoint(2, sf::Vector2f(w, P.L - P.b + h ));
    T_up.shape.setPoint(3, sf::Vector2f(0., P.L - P.b + h ));
    
    // up 
    BFace F_up;
    F_up.idx = N;
    F_up.iT_left = N-1;
    F_up.iT_right = N-1;
    F_up.xF = Vec::Zero(d);
    F_up.xF(0) = 0.5;
    F_up.xF(1) = P.b;
    F_up.xa = Vec::Zero(d);
    F_up.xb = Vec::Zero(d);
    F_up.xa(0) = 0.;
    F_up.xa(1) = P.b;
    F_up.xb(0) = w;
    F_up.xb(1) = P.b;
    F_up.vol = w;
    F_up.nF = Vec::Zero(d);
    F_up.nF(0) = 0;
    F_up.nF(1) = 1;
    F_up.pw = P.pw_top;
    F_up.sat = P.s_top;
    F_up.DIRICHLET_BC = true;
    
    // up left
    BFace F_up_l;
    F_up_l.idx = 2*N;
    F_up_l.iT_left = N-1;
    F_up_l.iT_right = 0;
    F_up_l.xF = Vec::Zero(d);
    F_up_l.xF(0)= 0.;
    F_up_l.xF(1)= P.b - h/2;
    F_up_l.xa = Vec::Zero(d);
    F_up_l.xb = Vec::Zero(d);
    F_up_l.xa(0) = 0.;
    F_up_l.xa(1) = P.b;
    F_up_l.xb(0) = 0.;
    F_up_l.xb(1) = P.b - h;
    F_up_l.vol = h;
    F_up_l.nF = Vec::Zero(d);
    F_up_l.nF(0) = -1;
    F_up_l.nF(1) = 0;
    F_up_l.NEUMANN_BC = true;
    
    BFace F_up_r;
    F_up_r.idx = 3*N;
    F_up_r.iT_left = N-1;
    F_up_r.iT_right = 0;
    F_up_r.xF = Vec::Zero(d);
    F_up_r.xF(0)= w;
    F_up_r.xF(1)= P.b-h/2;
    F_up_r.xa = Vec::Zero(d);
    F_up_r.xb = Vec::Zero(d);
    F_up_r.xa(0) = w;
    F_up_r.xa(1) = P.b;
    F_up_r.xb(0) = w;
    F_up_r.xb(1) = P.b - h;
    F_up_r.vol = h;
    F_up_r.nF = Vec::Zero(d);
    F_up_r.nF(0) = 1.;
    F_up_r.nF(1) = 0.;
    F_up_r.NEUMANN_BC = true;
    
    Mailles.push_back(T_bottom);    
    Mailles.push_back(T_up);
    BoundaryFaces.push_back(F_up);
    BoundaryFaces.push_back(F_up_l);
    BoundaryFaces.push_back(F_up_r);
    BoundaryFaces.push_back(F_bottom);
    BoundaryFaces.push_back(F_bottom_l);
    BoundaryFaces.push_back(F_bottom_r);
  } // boundary elements

  xmin = 0.;   xmax = w;
  ymin = 0.;   ymax = P.L;
  
  std::cout << "[Column2DMesh] N " << n  << std::endl;
  std::cout << "[Column2DMesh] width " << w  << std::endl;
  std::cout << "[Column2DMesh] xmin xmax " << xmin << " " << xmax << std::endl;
  std::cout << "[Column2DMesh] ymin ymax " << ymin << " " << ymax << std::endl;
  
  // sorting vectors
  std::sort(Mailles.begin(), Mailles.end(), compare<Element>);
  // not necessary to sort faces (?)
  //std::sort(FacesStd.begin(), FacesStd.end(), compare<Face>);
  //std::sort(BoundaryFaces.begin(), BoundaryFaces.end(), compare<BFace>);
  //std::sort(Interfaces.begin(), Interfaces.end(), compare<IFace>);
  
  ELEMENTS_ARE_SORTED = true;
}
