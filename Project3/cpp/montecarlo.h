#ifndef MONTECARLO_H
#define	MONTECARLO_H

#include <cmath>
#include <iostream>
#include <chrono>
using namespace std;
using namespace std::chrono;

void mc_bruteforce(double(*)(double *),int,double,double&,double&,double&);

void mc_improved(double(*)(double *),int,double&,double&,double&,double&,long);

#endif /* MONTECARLO_H */
