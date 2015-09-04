#include <iostream>
#include <Eigen/Dense>
#include "selfenergy.h"
#include "globals.h"
#include <list>
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <complex>
#include <mpi.h>


using namespace std;
using Eigen::VectorXcd;
using Eigen::VectorXd;

int main(int argc, char* argv[]){
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  if(myrank==0){
    cout << "starting program" << endl;
    cout << "numprocs" << numprocs << endl;
  }
  
  time_t tStart = time(NULL);
 
  kx = Q/2;
  ky = 1.94643;
  Nkpoints = 300;
  wLower = E(0,M_PI);
  wUpper = 0;
  Nf = 50;

  VectorXcd Epg = VectorXcd::Zero(2*Nf-1);
  double norm;
  double norm_private;
  VectorXcd one = VectorXcd::Ones(2*Nf-1);
  
  VectorXd wd = VectorXd::LinSpaced(2*Nf-1,wLower,wUpper);
  VectorXcd w = VectorXcd::Zero(2*Nf-1);
  for(int i=0; i<2*Nf-1; i++)
    w[i]=wd[i];


  while(true){
    
    if (myrank==0){
      cout << "Enter new D0_pg \n" << endl;
      cin >> D0_pg;
      D0_pg /= tval;
    }

    MPI_Bcast(&D0_pg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
  
    if (D0_pg == 0){
	  MPI_Finalize();
	  return 0;
    }
    

    MPI_Barrier(MPI_COMM_WORLD);

  if (1==1){ //pseudogap
    Epg = VectorXcd::Zero(2*Nf-1);
    norm = 0;
    norm_private = 0;
    

    if(myrank==0)
      cout << endl << "Calculating pseudogap selfenergies..." << endl;
  
    MPI_Barrier(MPI_COMM_WORLD);

    //calc_phonon_selfenergy(kx, ky, Eph_private, Ephm_private);
    calc_pseudogap_selfenergy(kx,ky,Epg,norm_private);
      
    //set up double[] for reduction
    complex<double> Epg_private_array[2*Nf-1];
    complex<double> Epg_array[2*Nf-1];
    for(int i=0;i<2*Nf-1;i++){
      Epg_private_array[i] = Epg[i];
    }

    MPI_Reduce(&Epg_private_array, &Epg_array, 2*Nf-1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&norm_private, &norm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(myrank==0)
      cout << "the norm is: " << norm << endl;

    if (myrank==0){   
      for(int i=0;i<2*Nf-1;i++){
	Epg[i] = Epg_array[i];
      }
      Epg /= norm;

      cout<<"Time taken "<< time(NULL)-tStart<<endl;

      ofstream myfile;
      myfile.open("Epg.dat");
      for (int i=0; i<2*Nf-1; i++){
	myfile << real(Epg[i]) << " " << imag(Epg[i]) << endl;
      }
      myfile.close();
      
      VectorXcd G = VectorXcd::Zero(2*Nf-1);
      //VectorXd A = VectorXd::Zero(2*Nf-1);
      double A[2*Nf-1]; 
      G = one.cwiseQuotient(w + (1j*gammae + E(kx,ky))*one - Epg);
      ///A = -imag(G);
      for (int i=0;i<2*Nf-1;i++){
      	 A[i] = -imag(G[i]);
	 cout << A[i] << endl;
      }

      double* biggest;
      double* first = A;
      double* last = first + 2*Nf-1;
      biggest = max_element(first,last);
      cout  << "freq at max " << w[distance(first,biggest)]*tval << endl;
      
      //double maxIndex;
      //G.maxCoeff(&maxIndex);
      
      //cout << maxIndex << endl;

      cout << endl <<  "Done With Program" << endl;
    }
  }

  }
  
  MPI_Finalize();
  return 0;
}
