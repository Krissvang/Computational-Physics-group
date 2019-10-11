#include <cmath>
#include <iostream>
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
  double int_mc, std_dev, time;
  mc_improved(&improved_MC,n,int_mc,std_dev,time);
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
