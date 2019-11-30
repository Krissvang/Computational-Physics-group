#ifndef VARMONTECARLO_H
#define	VARMONTECARLO_H

#include <armadillo>
using namespace std;
using namespace arma;

void var_mc(double&, double&, double&, int &, int, double, double, double,
            double (*)(mat &, double, double, double),
            double (*)(mat &, double, double, double));

mat init_pos();

double FindOptimal_h(double, double, double, int,
                     double (*)(mat &, double, double,double),
                     double (*)(mat &, double, double, double));

void solver(double&, double&, double&, int&, int, double, double, double,
            double (*)(mat &, double, double, double),
            double (*)(mat &, double, double, double), double);

double r_12(mat&);

double radi(int, mat&);

double r_squared(mat&);

double TrialWaveFunction1(mat&, double, double, double);

double TrialWaveFunction2(mat&, double, double, double);

double E1(mat&, double, double, double);

double E_repuls(mat&, double, double, double);

double E2(mat&, double, double, double);

double Kinetic_E(mat&, double, double, double, double, double, int,
                 double (*)(mat &r, double, double, double));

#endif /* VARMONTECARLO_H */
