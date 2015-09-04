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
 
  ifstream infile("pg.inp");
  read_params(infile);

  if (myrank==0){
    cout << endl << "##########  Params ########" << endl;
    cout << "D0_pg " << D0_pg*tval <<endl;
    cout << "gammacdw " << gammacdw << endl;
    cout << "gammae " << gammae*tval << endl;
    cout << "kx " << kx << endl;
    cout << "ky " << ky << endl;
    cout << "Nk " << Nkpoints << endl;
    cout << "Nf " << Nf << endl;
  }
  
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

      ofstream myfile;
      char buff [50];
      sprintf(buff,"Epg%1.2f.dat",D0_pg*tval);
      myfile.open(buff);
      //myfile.open("Epg.dat");
      for (int i=0; i<2*Nf-1; i++){
	myfile << real(Epg[i]) << " " << imag(Epg[i]) << endl;
      }
      myfile.close();

      cout << endl <<  "Done With Program" << endl;
    }
 
  
  MPI_Finalize();
  return 0;
}
