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

  ofstream myfile;

  if(myrank==0){
      myfile.open("EpgEk.dat");
  }
 
  int Ngrid = 60;
  for (int ikx=0; ikx<=Ngrid; ikx++){
    for(int iky=0; iky<=Ngrid; iky++){
 
      kx = ikx*M_PI/Ngrid;
      ky = iky*M_PI/Ngrid;
  
  if (1==1){ //pseudogap
    Nf = 1;
    wLower = E(kx,ky);
    wUpper = wLower;


    VectorXcd Epg = VectorXcd::Zero(2*Nf-1);
    double norm = 0;
    double norm_private = 0;
    
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

      for (int i=0; i<2*Nf-1; i++){
	myfile << real(Epg[i]) << " " << imag(Epg[i]) << " " << E(kx,ky) << endl;
      }
     
      cout << endl <<  "Done With Program" << endl;
    }
  } //if pseudogap

    } //iky loop
  } //ikx loop

  if (myrank==0)
     myfile.close();
  
  MPI_Finalize();
  return 0;
}
