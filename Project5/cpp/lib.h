#ifndef LIB
#define LIB

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

double r_12(mat &);
double radi(int, mat &);
double TrialWaveFunction1(mat &, double, double, double);
double TrialWaveFunction2(mat &, double, double, double);
double FindOptimal_h(double alpha, double, double, int,
                     double (*TrialFunction)(mat &, double, double, double));
double r_squared(mat &r);
double E1(mat &, double, double);
double E_repuls(mat &, double, double);
double E2(mat &, double, double,double);
double P1(mat &, double, double);
double P2(mat &, double, double,double);
mat init_pos();
double Kinetic_E(mat &, double, double, double, double, double, int,
                 double (*trialFunction)(mat &r, double, double, double));

#endif /* JACOBI_H */
