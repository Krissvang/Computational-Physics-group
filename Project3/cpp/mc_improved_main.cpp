#include <cmath>
#include <iostream>
#include <iomanip>
#include "lib.h"
#include "montecarlo.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

double improved_MC(double *);

int main(int argc, char const *argv[]) {
  int n;
  cout << "Read in the number of Monte-Carlo samples" << endl;
  cin >> n;
  n=pow(10,n);
  double int_mc, std_dev, time, sum_sigma;
  time_point<system_clock> time2;
  time2 = system_clock::now();
  duration<double> duration_in_seconds =duration<double>(time2.time_since_epoch());
  long t2= duration_in_seconds.count();
  mc_improved(&improved_MC,n,int_mc,std_dev,time,sum_sigma,t2);
  cout << "Standard deviation = "<< std_dev <<  " Integral = " << setprecision(10) << int_mc << " exact= " <<setprecision(10)<< 5*M_PI*M_PI/(16*16) << " Time = " << time << endl;
  return 0;
}

// this function defines the integrand to integrate
double  improved_MC(double *x)
{
   double cosb = cos(x[2])*cos(x[3])+sin(x[2])*sin(x[3])*cos(x[4]-x[5]);
   double deno = sqrt(x[0]*x[0]+x[1]*x[1]-2*x[0]*x[1]*cosb);
   double value = x[0]*x[0]*x[1]*x[1]*sin(x[2])*sin(x[3])/deno;
   return value;

} // end function for the integrand
