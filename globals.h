#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using Eigen::VectorXd; 
using namespace std;

#ifndef GLOBALS_H
#define GLOBALS_H

extern int myrank;
extern int numprocs;
extern double kx;
extern double ky;

extern double tval;
extern double t;
extern double t1;
extern double t2;
extern double t3;
extern double t4;
extern double tz;
extern double u;
extern double l_width;
extern int Nf;
extern double wLower;
extern double wUpper;
extern double Temperature;
extern double kT;
extern int Nkpoints;
extern int Nkz;
extern double gammacdw;
extern double Q;
extern double gammae;
extern double alpha;
extern double gamma_pg;
extern double T_star;
extern double D0_pg;

extern VectorXd kxs;
extern VectorXd kys;
extern int num_k_slices;
extern double kxmin;
extern double kxmax;
extern double kymin;
extern double kymax;


void read_item(int &x, ifstream &infile);
void read_item(double &x, ifstream &infile);
void read_params(ifstream &infile);


#endif
