#include <Eigen/Dense>
#include <vector>
#include <valarray>
#include "globals.h"
#include <complex>
using namespace std;
using Eigen::VectorXcd;
using Eigen::VectorXd;

#ifndef SELFENERGY_H
#define SELFENERGY_H

void calc_pseudogap_selfenergy(double kx, double ky, VectorXcd &Epg, double &norm);

double E(double kx, double ky);

double PseudoGap();

double PseudoGap(double kx,double ky);

#endif
