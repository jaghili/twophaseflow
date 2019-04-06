#include "defs.h"
#include "FVCA5Mesh.h"
#include <Eigen/Dense>
#include "Element.h"
#include "Face.h"
#include "BFace.h"
#include <fstream>
#include <SFML/Graphics.hpp>

FVCA5Mesh::FVCA5Mesh(const char* file, R scale) {
  std::fstream datafile(file, std::ios::in);
  std::string marker;

  // fixed dim
  d = 2;
  Nrt = 2;
  
  // vertices
  datafile >> marker;
  while(marker.compare("vertices")) {
    datafile >> marker;
    assert(!datafile.eof());
  }
  uint n_vertices;
  datafile >> n_vertices;
  std::cout << "[FVCA5 Mesh] Reading mesh " << file << std::endl;
  std::cout << "[FVCA5 Mesh] Scale factor =" << scale << "\n";  
  std::cout << "[FVCA5 Mesh] Reading " << n_vertices << " vertices...\n";
  Mat points(d, n_vertices);
  for (uint i = 0; i < n_vertices; i++)
    {
      Vec v = Vec::Zero(d);
      datafile >> v(0) >> v(1);
      //std::cout << "Point no. " << i << " \n"<< v << std::endl;
      points.col(i) = scale * v;
    }
  xmin = points.row(0).minCoeff();   xmax = points.row(0).maxCoeff();
  ymin = points.row(1).minCoeff();   ymax = points.row(1).maxCoeff(); 
  std::cout << "[FVCA5 Mesh] Dimensions x=" << xmin << "x" << xmax << "\t y=" << ymin << "x" << ymax << "\n" ;
  //
  N=0;
  
  // triangles elements
  uint n_triangles = 0;
  datafile >> marker;
  while(marker.compare("triangles")) {
    datafile >> marker;
    assert(!datafile.eof());
  }
  datafile >> n_triangles;
  N += n_triangles;
  std::cout << "[FVCA5 Mesh] Reading " << n_triangles << " triangles...\n";  
  Mailles.reserve(N);
  for (uint i = 0; i < n_triangles; i++) {
    uint iv1,iv2,iv3;
    datafile >> iv1 >> iv2 >> iv3;

    Element T;
    
    T.sommets = Eigen::MatrixXd::Zero(d,4);
    T.sommets.col(0) = points.col(iv1-1);
    T.sommets.col(1) = points.col(iv2-1);
    T.sommets.col(2) = points.col(iv3-1);
    T.sommets.col(3) = points.col(iv1-1);
        
    T.idx = i;
    T.vol = 0.5 * std::abs(
			   T.sommets.block<2,2>(0,0).determinant()
			   + T.sommets.block<2,2>(0,1).determinant()
			   + T.sommets.block<2,2>(0,2).determinant()
			   );
    T.xT = Vec::Zero(d);
    T.xT = (1./3.) * (T.sommets.col(0) + T.sommets.col(1) + T.sommets.col(2));
    T.rt = 0;
    T.shape.setPointCount(3);
    for (int j = 0; j < 3; j++) {
      T.shape.setPoint(j, sf::Vector2f(T.sommets(0,j), ymax - T.sommets(1,j)));
    }

    Mailles.push_back(T);
  }

  // quadrangle elements
  uint n_quadrangles = 0;
  datafile >> marker;
  while(marker.compare("quadrangles")) {
    datafile >> marker;
    assert(!datafile.eof());
  }
  datafile >> n_quadrangles;
  N += n_quadrangles;
  std::cout << "[FVCA5 Mesh] Reading " << n_quadrangles << " quadrangles...\n";
  // Adding elements
  Mailles.reserve(N);
  for (uint i = 0; i < n_quadrangles; i++) {
    uint iv1,iv2,iv3,iv4;
    datafile >> iv1 >> iv2 >> iv3 >> iv4;

    Element T;
    
    T.sommets = Eigen::MatrixXd::Zero(d,5);
    T.sommets.col(0) = points.col(iv1-1);
    T.sommets.col(1) = points.col(iv2-1);
    T.sommets.col(2) = points.col(iv3-1);
    T.sommets.col(3) = points.col(iv4-1);
    T.sommets.col(4) = points.col(iv1-1);    

    T.idx = i;
    T.vol = 0.5 * std::abs(
			   T.sommets.block<2,2>(0,0).determinant()
			   + T.sommets.block<2,2>(0,1).determinant()
			   + T.sommets.block<2,2>(0,2).determinant()
			   + T.sommets.block<2,2>(0,3).determinant()
			   );
    T.xT = Vec::Zero(d);
    T.xT = 0.25 * (T.sommets.col(0) + T.sommets.col(1) + T.sommets.col(2) + T.sommets.col(3));
    T.rt = 0;
    T.shape.setPointCount(4);
    for (int j = 0; j < 4; j++) {
      T.shape.setPoint(j, sf::Vector2f(T.sommets(0,j), ymax - T.sommets(1,j)));
    }
    //std::cout << "Element no." << i << " center=" << T.xT << std::endl;
    //std::cout << "Adding element no. " << T.idx << " v1=" << iv1 << " v2=" << iv2 << " v3=" << iv3 << " v4=" << iv4 << std::endl; 
    Mailles.push_back(T);
  }
  
  // boundary faces
  datafile >> marker;
  while(marker.compare("boundary")) {
    datafile >> marker;
    assert(!datafile.eof());
  }
  datafile >> Nfb;
  std::cout << "[FVCA5 Mesh] Reading " << Nfb << " boundary faces...\n";
  
  // all faces
  datafile >> marker;
  while(marker.compare("edges")) {
    datafile >> marker;
    assert(!datafile.eof());
  }
  uint n_faces;
  datafile >> n_faces;
  FacesStd.reserve(n_faces - Nfb);
  BoundaryFaces.reserve(Nfb);
  std::cout << "[FVCA5 Mesh] Reading " << n_faces << " faces...\n";
  uint intF = 0;
  uint bdF = 0;
  for (uint iF = 0; iF < n_faces; iF++) {
    uint iv1, iv2;
    uint iTl, iTr;
    datafile >> iv1 >> iv2 >> iTl >> iTr;
    //std::cout << "Face no. " << iF << " vertices=" << iv1 << " " << iv2 << "\n";
    Vec _xF = 0.5 * (points.col(iv1-1) + points.col(iv2-1));
    R   _vol =  (points.col(iv1-1) - points.col(iv2-1)).norm();
    assert(iv1-1 < n_vertices and iv2-1 < n_vertices);
    
    if (iTr == 0 and iTl > 0) { // boundary face
      bdF++;
      BFace F;
      F.idx = iF;
      F.iT_left = iTl-1;
      F.iT_right = iTl-1;
      F.vol = _vol;
      F.xF = _xF;
      F.xa = points.col(iv1-1);
      F.xb = points.col(iv2-1);
      F.nF = (_xF - Mailles[iTl-1].xT) / (_xF - Mailles[iTl-1].xT).norm();
      // addind the face
      BoundaryFaces.push_back(F);
    }

    else if (iTr > 0 and iTl > 0) { // interior face
      intF++;
      Face F;
      F.idx = iF;
      F.iT_left = iTl-1;
      F.iT_right = iTr-1;
      F.vol = _vol;
      F.xF = _xF;
      F.xa = points.col(iv1-1);
      F.xb = points.col(iv2-1);
      F.nF = (Mailles[iTr-1].xT - Mailles[iTl-1].xT) / (Mailles[iTr-1].xT - Mailles[iTl-1].xT).norm(); // vecteur normal de left vers right
      // adding the face 
      FacesStd.push_back(F);
    }
  }
  Nfb =  bdF;
  Nf = intF;

  // rocktype
  datafile >> marker;
  if (!marker.compare("rocktype")) {
    std::cout << "[FVCA5 Mesh] Setting rocktypes..." << std::endl;
    uint no_iT = 0;
    while(!datafile.eof()) {
      datafile >> Mailles[no_iT].rt;
      no_iT++;
    }
  }
    
  // sorting Mailles
  std::sort(Mailles.begin(), Mailles.end(), compare<Element>);
  ELEMENTS_ARE_SORTED = true;

  //
  datafile.close();
}
