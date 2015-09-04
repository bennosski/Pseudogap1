#include <Eigen/Dense>
#include "selfenergy.h"
#include "globals.h"
#include <math.h>
#include <complex>
#include <vector>
#include <iostream>
#include <omp.h>
#include <mpi.h>
using namespace std;
using Eigen::VectorXcd;
using Eigen::VectorXd;


void calc_pseudogap_selfenergy(double kx,double ky, VectorXcd &Epg, double &norm){
   MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
   int namelen;
   char processor_name[MPI_MAX_PROCESSOR_NAME];
   MPI_Get_processor_name(processor_name, &namelen);

   cout << "hello rank " << myrank << " " << processor_name << endl;
   
   VectorXcd one = VectorXcd::Ones(2*Nf-1);
   VectorXd wd = VectorXd::LinSpaced(2*Nf-1,wLower,wUpper);
   VectorXcd w = VectorXcd::Zero(2*Nf-1);

   for(int i=0; i<2*Nf-1; i++)
     w[i]=wd[i];

   int Nk = Nkpoints;
   double totalN = pow(2*Nk,2)*(2*Nkz);
   if (Nkz==1)
     totalN = pow(2*Nk,2);
   double pgap = PseudoGap(kx,ky);
   complex<double> i1(0.0,1.0);

   double qx,qy,qz,px,py,pofq,Energy;

      for (int q1=-Nk+myrank; q1<Nk; q1+=numprocs){ //BZ momentum loop      
         qx = M_PI*q1/Nk;
	 for(int q2=-Nk; q2<Nk; q2++){
	    qy = M_PI*q2/Nk;
	    
	    if (Nkz==1){
	      px = kx - qx;
              py = ky - qy;

	      Energy = E(px,py);

              pofq = 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx-Q,2)+pow(qy,2)) + 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx+Q,2)+pow(qy,2)) + 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx,2)+pow(qy-Q,2)) + 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx,2)+pow(qy+Q,2));
	    
	      //pofq = 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx-Q,2)+pow(qy,2)) + 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx+Q,2)+pow(qy,2));


	      Epg += pgap*pofq*one.cwiseQuotient(w-one*E(px,py)+i1*gammae*one)/totalN;
	      norm += pofq/totalN;
	    }
	    else{
	    for(int q3=-Nkz; q3<Nkz; q3++){
	      qz = M_PI*q3/Nkz;	

	      px = kx - qx;
              py = ky - qy;

              pofq = 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx-Q,2)+pow(qy,2)) + 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx+Q,2)+pow(qy,2)) + 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx,2)+pow(qy-Q,2)) + 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx,2)+pow(qy+Q,2));
	    
	      //pofq = 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx-Q,2)+pow(qy,2)) + 1.0/M_PI*gammacdw/(pow(gammacdw,2)+pow(qx+Q,2)+pow(qy,2));

	      Epg += pgap*pofq*one.cwiseQuotient(w-one*(E(px,py)-2*tz*cos(qz))+i1*gammae*one)/totalN;
	      norm += pofq/totalN;
	    }
	    }

         }
      } // end BZ loop

}

double E(double kx, double ky){
  //norman bandstructure
  return u+0.5*t*(cos(kx)+cos(ky))+t1*cos(kx)*cos(ky)+0.5*t2*(cos(2*kx)+cos(2*ky))+0.5*t3*(cos(2*kx)*cos(ky)+cos(kx)*cos(2*ky))+t4*cos(2*kx)*cos(2*ky);
 
  //Moritz bandstructure
  //return -2.0*t*(cos(kx)+cos(ky))-4*t1*cos(kx)*cos(ky)-2*t2*(cos(2*kx)+cos(2*ky))-4*t3*(cos(2*kx)*cos(ky)+cos(kx)*cos(2*ky))+u;
    

  //return -2*t*(cos(kx) + cos(ky)) - 4*t1*cos(kx)*cos(ky) - u;
}

double PseudoGap(){
  return 0;
}

double PseudoGap(double kx,double ky){
  return D0_pg;
  //return PseudoGap(); //*pow((cos(kx)-cos(ky))/2.0,3);
}

