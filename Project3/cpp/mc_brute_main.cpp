#include <cmath>
#include <iostream>
#include "lib.h"
#include "montecarlo.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

double brute_force_MC(double *);

int main(int argc, char const *argv[]) {
  int n;
  double R;
  cout << "Read in the number of Monte-Carlo samples" << endl;
  cin >> n;
  cout << "Read in R (max. absolute value of x_i)" << endl;
  cin >> R;
  double int_mc, std_dev, time;
  mc_bruteforce(&brute_force_MC,n,R,int_mc,std_dev,time);
  return 0;
}

// this function defines the integrand to integrate
double  brute_force_MC(double *x)
{
   double alpha = 2.;
// evaluate the different terms of the exponential
   double exp1=-2*alpha*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
   double exp2=-2*alpha*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);;
   double deno = sqrt(pow((x[0]-x[3]),2)+pow((x[1]-x[4]),2)+pow((x[2]-x[5]),2));
   return exp(exp1+exp2)/deno;
} // end function for the integrand
