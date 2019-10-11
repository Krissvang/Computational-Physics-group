
#include <cmath>
#include <iostream>
#include "lib.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

double brute_force_MC(double *);
//     Main function begins here
int main()
{
     int n;
     double x[6], y, fx;
     cout << "Read in the number of Monte-Carlo samples" << endl;
     cin >> n;
     time_point<high_resolution_clock> start, end;
     start = high_resolution_clock::now();
     double int_mc = 0.;  double variance = 0.;
     double sum_sigma= 0. ; long idum=-1 ;
     double length=10.; // we fix the max size of the box to L=20
     double jacobi_det=pow((2*length),6);

//   evaluate the integral with importance sampling
     for ( int i = 1;  i <= n; i++){
//   x[] contains the random numbers for all dimensions
       for (int j = 0; j< 6; j++) {
           x[j]=-length+2*length*ran0(&idum);
       }
       fx=brute_force_MC(x);
       int_mc += fx;
       sum_sigma += fx*fx;
     }
     int_mc = int_mc/((double) n );
     sum_sigma = sum_sigma/((double) n );
     variance=sum_sigma-int_mc*int_mc;
     end = high_resolution_clock::now();
     duration<double> elapsed = end-start;
     double time = elapsed.count();
//   final output
     cout << "Standard deviation=" << jacobi_det*sqrt(variance/n) <<" Integral= " << jacobi_det*int_mc <<" Exact= " << 5*M_PI*M_PI/(16*16) << " Time = "<< time << endl;
     return 0;
}  // end of main program

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
