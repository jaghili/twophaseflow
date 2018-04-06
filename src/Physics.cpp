#include "Physics.h"
#include <fstream>

void Physics::writetofile() {
  std::ofstream pcsat("./output/physics_"+pcname+"_sat.dat");
  std::ofstream pcsati("./output/physics_"+pcname+"_tau.dat");
  std::ofstream pdata("./output/physics.tmp");
  
  // writing physical data
  pdata << "Physics \n"
	<< "a\t\t\t" << a << "\n"
	<< "b\t\t\t" << b << "\n"
	<< "L\t\t\t" << L << "\n"
	<< "rhow\t\t\t" << rhow << "\n"
	<< "rhoo\t\t\t" << rhoo << "\n"
	<< "grav\t\t\t" << grav << "\n"
	<< "sor\t\t\t" << sor << "\n"
	<< "swi\t\t\t" << swi << "\n"
	<< "muo\t\t\t" << muo << "\n"
	<< "muw\t\t\t" << muw << "\n"
	<< "pw\t\t\t" << pw << "\n"
	<< "phi\t\t\t" << phi << "\n"
	<< "Nrt\t\t\t" << Nrt << "\n"
	<< "tau1\t\t\t"  << tau1 << "\n"
    	<< "tau2\t\t\t"  << tau2 << "\n"
    	<< "tau3\t\t\t"  << tau3 << "\n"
    	<< "b1\t\t\t"  << b1 << "\n"
    	<< "b2\t\t\t"  << b2 << "\n"
	<< "perm\t\t\t" << perm[0] << "\t" << perm[1] << "\n\n";

  pdata << "Boundary Conditions\n"
	<< "s_bottom\t\t\t" << s_bottom << "\n"
	<< "s_top\t\t\t" << s_top << "\n"
	<< "pw_inj\t\t\t" << pw_inj << "\n"
	<< "pw_bottom\t\t\t" << pw_bottom << "\n"
	<< "pw_top\t\t\t" << pw_top << "\n";
    
  // writing pc/sat laws
  // faces
  pcsat << "sat\tpc0\tdpc0\tpc1\tdpc1" << "\n";
  for (R sat = 0.; sat <= 1.; sat += 1e-2) {
    pcsat << sat
	  << "\t" << pc[0][0](sat)
	  << "\t" << dpc[0][0](sat)
	  << "\t" << pc[1][1](sat)
	  << "\t" << dpc[1][1](sat) << "\n";
  }
  {
    R sat = 1.;
    pcsat << sat
	  << "\t" << pc[0][0](sat)
	  << "\t" << dpc[0][0](sat)
	  << "\t" << pc[1][1](sat)
	  << "\t" << dpc[1][1](sat) << "\n";
  }
  
  // intefaces
  pcsati << "tau\tpc01\tdpc01\tsat01\tdsat01\tsat10\tdsat10" << "\n";
  for (R tau = taumin; tau < taumax; tau += 1e-2) {
    pcsati << tau
	   << "\t" << pc[0][1](tau)
	   << "\t" << dpc[0][1](tau)
	   << "\t" << sat[0][1](tau)
	   << "\t" << dsat[0][1](tau)
	   << "\t" << sat[1][0](tau)
	   << "\t" << dsat[1][0](tau) << "\n";
  }
  {
    R tau = taumax;
    pcsati << tau
	   << "\t" << pc[0][1](tau)
	   << "\t" << dpc[0][1](tau)
	   << "\t" << sat[0][1](tau)
	   << "\t" << dsat[0][1](tau)
	   << "\t" << sat[1][0](tau)
	   << "\t" << dsat[1][0](tau) << "\n";
  }
  pcsati.close();
  pcsat.close();
  pdata.close();
}
