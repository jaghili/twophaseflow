#define EIGEN_STACK_ALLOCATION_LIMIT 0

#include <iostream> // IO
#include <iomanip> // std::scientific stuff
#include <fstream> // file IO
#include <functional> // to define lambda functions
#include <vector>
#include <random>

// for live plotting 
#include <SFML/Window.hpp>
#include <SFML/Graphics.hpp>

#include <Eigen/Core>
#include <Eigen/SparseLU> // Sparse LU (not superLU!)

#define VIENNACL_WITH_EIGEN 1

//#include <Eigen/UmfPackSupport>
//#include <Eigen/SuperLUSupport>
//#include <Eigen/PardisoSupport>
// Iterative solvers
//#include <Eigen/IterativeSolvers> // GMRES in unsupported
//#include <Eigen/IterativeLinearSolvers> // BiCGSTAB


// GPU support via viennacl
//#include "viennacl/scalar.hpp"
//#include "viennacl/vector.hpp"
//#include "viennacl/compressed_matrix.hpp"
//#include "viennacl/coordinate_matrix.hpp"
//#include "viennacl/linalg/prod.hpp"
//#include "viennacl/linalg/gmres.hpp"
//#include "viennacl/linalg/ilu.hpp"


#include "GetPot" // library to parse arguments 

#include "Mesh.h" 
#include "Physics.h"
#include "Corey.h"
#include "Constants.h"
#include "CoreyConstant.h"
#include "Column.h"
#include "FVCA5Mesh.h"
#include "Data.h"
#include "Residual.h"
#include "Jacobian.h"
#include "DrawSolution.h"
#include "InterfaceSolver.h"
#include "BasicFunctions.h"


// ---- I/O -----------
std::ofstream nwt_file("./output/newton2d_convergences.dat");     // newton iterations for the full problem

// "extern" tells the compiler that computeF() exists and is defined outside of this file (in Flux.cpp)
extern void computeF(Mesh* Mh, Dat & Xold, Physics & P);

// ---------- Newton settings ----------------
size_t nwt_iter               =0;
size_t nwt_global_iters       =0;
float nwt_eps                 =1e-11;
float nwt_dxmax               =1e-7;
uint nwt_max_iter             = 100;
uint nwt_is_fails             = 0;
uint nwt_fails                = 0;
float nwt_dSobj               = 0.25;
R dxmax                       = 1e100;
R nwt_offset                  = 0.;
R nwt_offset_iface            = 0.;

// ---- Time settings (temps en secondes)  -----------------------
uint nt                       = 0;
R Temps                       = 0;
R UneHeure                    = 3600.; 
R UnJour                      = 24 * UneHeure;
R UnAn                        = 360 * UnJour;
R Tmax                        = 5000 * UnAn; 
R dtinit                      = UnJour; // 1 jours (s)
R dtmax                       = UnAn*100; // 10 ans (s)
R dt                          = dtinit;

bool use_is                   = true;
R vvol                        = 0.0001;


// ---
sf::Color c_white(220,220,220);

int main(int argc, char * argv[]) {

  // get arguments
  GetPot options(argc, argv);

  // Choose capillary pressure law
  //CoreyConstant P; // P contains all the physics constants (with Constants laws)
  //Constants P; // P contains all the physics constants (with Constants laws)
  Corey P; // P contains all the physics constants (with Corey laws)
  
  std::cout << "[Physics] Pc laws\t" << P.pcname << "\n";
  P.writetofile();
  //exit(EXIT_SUCCESS);
  
  // FVCA5 Mesh  
  std::string mesh_file = "meshes/fracturemaille.msh";
  R scale = 1.;
  if(options.search(2, "--mesh", "-m")) { mesh_file = options.next(""); }
  if(options.search(2, "--scale", "-s")) { scale = stod(options.next("")); }
  
  FVCA5Mesh Mh(mesh_file.c_str(), scale);
  //Column Mh(101, 100., P); // Build a 2D Column Mesh using Physics  
  sf::View view(sf::FloatRect(Mh.xmin, Mh.ymin, Mh.xmax, Mh.ymax));
  
  // Plot : Render window 
  sf::RenderWindow window(sf::VideoMode(800,800), "Plot");
  
  //window.setFramerateLimit(30); // FPS limiter
  window.setView(view); // focus the camera on the (0,100)x(0,100) region of the window

#ifdef RECORDMOVIE
  sf::Texture texture;
  sf::Vector2u windowSize = window.getSize();
  texture.create(windowSize.x, windowSize.y);
#endif
  
  // Rocktype initialization
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0,1);
  for (auto & M : Mh.Mailles) {
    R r = dis(gen);
    //M.rt = ( (r < 0.2 and M.xT(1) > 10) ? 1 : 0); // 1 : barrière , 0: drain
    //M.rt = ((M.xT(1) > 10 and M.xT(1) < 20) or (M.xT(1) > 25 and M.xT(1) < 30)? 1 : 0 );
    //M.rt = (M.xT(1) > 10 and M.xT(1) < 15 ? 1 : 0 );
  }
  
  // Compute interfaces according to rocktype distribution
  Mh.computeInterfaces(); // interfaces

  
  // Boundary conditions 
  uint DirBC = 0;
  uint NeuBC = 0;
  for (auto & F : Mh.BoundaryFaces) {    
    if (std::abs(F.xF(1)) < 1. and std::abs(F.xF(0) - 25.) < F.vol/2) {
    //if (std::abs(F.xF(1)) < 1. and std::abs(F.xF(0) - 50.) < F.vol) {
    //if (std::abs(F.xF(1)) < 1.) { // DOWN
      F.setDirichletBC(P.s_bottom, P.pw_bottom);
      DirBC++;
    }
    else if (std::abs(F.xF(1) - P.L) < 1.) { // UP
      F.setDirichletBC(P.s_top, P.pw_top);
      DirBC++;
    } // if F is on the top
    else {
      F.setNeumannBC(0., 0., 0., 0.);
      NeuBC++;
    } // elsewhere 
  }
  
  // Build solution vector
  Dat Xold(&Mh);
  
  // Initial conditions
  Xold.S = Vec::Zero(Mh.N);
  Xold.P = Vec::Zero(Mh.N);

  // Recap
  std::cout << "\n";
  std::string prefix = "[Settings] ";
  std::cout << prefix << "no. elements \t"<< Mh.N << std::endl;
  std::cout << prefix << "no. interfaces \t"<< Mh.Ni << std::endl;
  std::cout << prefix << "no. faces \t"<< Mh.Nf << std::endl;
  std::cout << prefix << "no. boundary faces \t"<< Mh.Nfb << std::endl;
  std::cout << prefix << "no. Dirichlet BC faces \t" << DirBC << std::endl;
  std::cout << prefix << "no. Neumann BC faces \t" << NeuBC << std::endl;
  std::cout << prefix << "no. dofs \t"<< Xold.Ndof << std::endl;
  std::cout << prefix << "no. interfaces \t"<< Mh.Interfaces.size() << "\n" << std::endl;

  R VolSat[P.Nrt] = { 0., 0. };
  // saturation initiale
  /*
  for (auto & M : Mh.Mailles) {
    std::cout << "vol maille " << M.vol << std::endl;
    //Xold.S(M.idx) = (M.xT(1) < 50.? 0. : 1.);
    //Xold.S(M.idx) = (M.xT(1) < 40. or M.xT(1) > 60 ? 1. : 0.5);
  }
  */
  // Pression initiale loi hydrostatique
  for (int i = 0 ; i < Xold.N; i++) {
    Xold.P(i) = P.pw + (P.L - Mh.Mailles[i].xT(1) ) * P.grav * P.rhow;
  }
  
  // Show Xold
  //Xold.display();

  // Show mesh
  //Mh.display();
  
  // Initialize solver(s) 
  Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int> > solver; // Eigen's Sparse LU
  //Eigen::SuperLU<SpMat> solver;                                 // SuperLU wrapper
  //Eigen::PardisoLU<SpMat> solver;                             // Intel's LU solver (super fast!)
  //Eigen::UmfPackLU<SpMat> solver;
  //Eigen::GMRES<SpMat, Eigen::IncompleteLUT<R> > solver; 
  /*
    Eigen::BiCGSTAB<SpMat, Eigen::IncompleteLUT<R> > solver;
    solver.setTolerance(1e-10);
    solver.preconditioner().setDroptol(1e-14);
    solver.preconditioner().setFillfactor(100);
  */
  
  // Init Solution container
  Dat X(Xold);
  std::string satfile = "./output/volumesat.dat";
  std::ofstream writesat(satfile);


  
  // ---------------------- Time Loop -------------------------------------- //
  std::cout << std::endl;
  while ( Temps < Tmax ) {

    // draw
    window.clear(c_white);
    DrawSolution(&Mh, P, X, window);
    window.display();


    // courbes d'huile/eau
    for (auto & M : Mh.Mailles) { VolSat[M.rt] += M.vol * X.S(M.idx) * P.phi[M.rt]; }
    writesat << Temps/UnJour << " " << VolSat[0] << " " << VolSat[1] << std::endl;
  
    // Compute Interface Values with IS
    if (use_is) {
      vvol = 0.;
      for (auto & F : Mh.Interfaces) {
	Vec Res = InterfaceSolver(&Mh, X, F, P, 1e-14, 1e-2, 1e-13);
      }
    }

    
    // snapshot to jpg
#ifdef RECORDMOVIE
    prefix = "[CAPTURE] ";
    std::string capturefile = "./output/movie/" + P.pcname + "_" + std::to_string(nt) + ".bmp";
    texture.update(window);
    texture.copyToImage().saveToFile(capturefile);
    std::cout << prefix << "screenshot in " << capturefile << std::endl;
#endif

    // Time update
    nt++;
    Temps += dt;

    // print on screen
    prefix = "[Time] ";
    std::cout << prefix << "n=" << nt << std::setprecision(15)
	      <<"\t t=" << Temps/UnJour
	      << "\tdt=" << dt/UnJour << "\n";


    // ---------------------- Newton Loop ------------------------------------ //
    prefix = "[Newton2D] ";
    R nwt_error = 1e300;
    nwt_iter=0;
    X = Xold;

    // GPU : prepare data 
    //viennacl::vector<double> vcl_B(2*X.Ndof);
    //viennacl::compressed_matrix<double> vcl_J(2*X.Ndof, 2*X.Ndof);
    
    std::cout << prefix
	      << std::setw(5) << std::scientific << "globaliters"
	      << std::setw(30) << std::scientific << "err"
	      << std::endl;

    
    while (nwt_error > nwt_eps and nwt_iter < nwt_max_iter) {
      // newton iter
      nwt_iter++;
      nwt_global_iters++;

      // Incrementation Newton      
      // compute flux
      computeF(&Mh, X, P); // 

      // compute residual and jacobian      
      auto Res = computeResidual(&Mh, X, Xold, P); // compute residual
      auto dRes = computeJacobian(&Mh, X, Xold, P); // compute jacobian of residual
      
      //viennacl::copy(dRes, vcl_J);
      //viennacl::linalg::ilu0_precond<viennacl::compressed_matrix<R>> vcl_ilu(dRes, viennacl::linalg::ilu0_tag()); // compute ilu0 prec with CPU
      
      nwt_error = Res.lpNorm<2>();
      
      // print on screen newton iterations
      std::cout << prefix
		<< std::setw(5) << std::scientific << nwt_global_iters
		<< std::setw(30) << std::scientific << nwt_error
		<< std::endl;
      
      // write newton stats 
      nwt_file << nwt_offset + nwt_iter << "\t" << nwt_error << std::endl;
      
      // solve
      Vec B = Vec::Zero(2*X.Ndof); // Rhs      
      for (uint i = 0; i < Res.rows(); i++) {
	B(2*i)   = -Res(i,0);
	B(2*i+1) = -Res(i,1);
      }
      // CPU factorization
      solver.compute(dRes); // analyze and factorize dRes

      // GPU factorize
      //viennacl::copy(B, vcl_B);

      /*
      if(solver.info() != Eigen::Success) {	
	// decomposition fails
	std::cout << prefix << "\033[7;33m Solver failed! Exporting matrix and RHS. Terminating...\033[0m\n";
	std::ofstream jac("jac");
	std::ofstream rhs("rhs");

	jac << dRes;
	rhs << B;
	
	rhs.close();
	jac.close();
	return 1;
      }
      */

      // compute increment
      Vec dU = Vec::Zero(2*X.Ndof);
      dU = solver.solve(B); // direct solver on CPU
      //viennacl::linalg::gmres_tag my_gmres_tag(1e-10, 50, 30); // up to 100 iterations, restart after 20 iterations
      //viennacl::vector<R> vcl_dU = viennacl::linalg::solve(vcl_J, vcl_B, my_gmres_tag, vcl_ilu); // on GPU
      //dU = viennacl::linalg::solve(dRes, B, my_gmres_tag); // on GPU
      //vcl_B = viennacl::linalg::solve(vcl_J, vcl_B, viennacl::linalg::bicgstab_tag());
      //std::cout << "\t gmres iters " << my_gmres_tag.iters() << "\t err " << my_gmres_tag.error() << std::endl;
      //viennacl::copy(vcl_dU, dU);
      
      R dSmax = 0;
      R dPmax = 0;
      for (uint i = 0; i < X.Ndof; i++) {
	dPmax = std::max(dPmax, std::abs(dU(2*i)));
	dSmax = std::max(dSmax, std::abs(dU(2*i+1)));
      }
      R alpha = std::min(1., nwt_dSobj/dSmax);
      dxmax = dSmax + dPmax * 1e-6;

      const Dat Xold = X;
      Vec vX2 = X.get() + alpha * dU;
      X.set(vX2);       // update X

      // Project variables
      P.project(X, Xold, P.tau1, P.tau2, P.tau3);
    }  // end newton loop

    
    // if newton has converged
    if (nwt_iter < nwt_max_iter) {
      //std::cout << "Newton \033[1;32m OK\033[0m \titers=" << nwt_iter << std::endl; 
      dt = (1.3*dt>dtmax ? dtmax : 1.3*dt); // control dt increase
      dt = (Temps+dt>Tmax ? Tmax-Temps : dt); // 
    } else {
      nwt_fails++;
      std::cout << prefix << "\033[7;33mFAILED\033[0m (too much iters) \titer=" << nwt_iter
		<< "\tRestarting..."<< std::endl;
      //exit(EXIT_FAILURE);
      Temps -= dt;
      nt -= 1;
      dt = 0.5*dt;
    }
    nwt_offset += 3;
    nwt_file << std::endl << std::endl;
    
    // update 
    Xold=X;

    std::cout << "\n" << LONGLINE << std::endl;
       
  } // end Time loop
  writesat.close();
  
  // writing solutions
  std::cout << std::endl;
  std::cout << LONGLINE << std::endl;
  std::cout << "Writing solution....\n";
  std::ofstream sol_file("./output/finalsol.tmp");        // store solutions
  for (uint i=0; i < Mh.N ; i++) {
    sol_file << std::setprecision(16) << std::scientific
	     << X.S(i) << "\t" // saturation
	     << X.P(i) << "\t" // wat pressure
	     << (i<Mh.Ni?X.Pi(i):0.) << "\t" // Pressure at interface
	     << (i<Mh.Ni?X.Tau(i):0.) 
	     << std::endl; 
  }
  sol_file.close();
  
  // Recap
  std::cout << LONGLINE << std::endl;
  std::cout << "Newton fails: \t" << nwt_fails << std::endl;
  // close window
  window.close();

  // close files
  nwt_file.close();

  // end
  std::cout << "Done.\n";
  return 0;
}


// Local Variables:
// compile-command: "make -C ../"
// End:
