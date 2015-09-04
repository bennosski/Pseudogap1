#include "globals.h"
#include <Eigen/Dense>
#include <math.h>
#include <fstream>
#include <string>
using namespace std;
using Eigen::VectorXd;

int myrank,numprocs;

double kx;
double ky;

//////norman bandstructure
double tval = 595.1/2;
double u = 120.0/tval;
double t = -595.1/tval;
double t1 = 163.6/tval;
double t2 = -51.9/tval;
double t3 = -111.7/tval;
double t4 = 51.0/tval;
double tz = 10/tval;


//////Moritz bandstructure
/*
double tval = 220.0;
double t = 220.0/tval;
double t1 = -34.315/tval;
double t2 = 35.977/tval;
double t3 = -7.1637/tval; 
double u = -169.0/tval;
double tz = 10/tval;
*/

int Nkpoints; 
int Nkz;
int Nf;
double D0_pg;
double wUpper;
double wLower;
double gammacdw;
double gammae;
//for 11.3 degree nesting
//double Q = 1.6828;
//for 18 degree nesting
//double Q = 1.43467;
double Q;

double l_width = 0.017;
double Temperature = 90;
double kT = .086173 * Temperature / tval; //units of t
double alpha = 1.0;
double gamma_pg = 0.0;
double T_star = 200.0;

int num_k_slices;
VectorXd kxs;
VectorXd kys;
double kxmin; double kxmax;
double kymin; double kymax;


void read_item(int &x, ifstream &infile){
  string nothing = "";
  infile >> nothing;
  infile >> x;
}

void read_item(double &x, ifstream &infile){
  string nothing = "";
  infile >> nothing;
  infile >> x;
}

void read_params(ifstream &infile){
  read_item(D0_pg, infile);  D0_pg /= tval;
  read_item(Nkz, infile);
  read_item(Nkpoints, infile);
  read_item(Nf, infile);
  read_item(kx, infile);
  read_item(ky, infile);
  read_item(gammacdw, infile);
  read_item(gammae, infile);  gammae /= tval;
  read_item(Q, infile);
  read_item(wLower, infile);  wLower /= tval;
  read_item(wUpper, infile);  wUpper /= tval;
  read_item(num_k_slices,infile);
  read_item(kxmin, infile);
  read_item(kxmax, infile);
  read_item(kymin, infile);
  read_item(kymax, infile);
  
  kxs = VectorXd::LinSpaced(num_k_slices,kxmin, kxmax);
  kys = VectorXd::LinSpaced(num_k_slices,kymin, kymax);
}
