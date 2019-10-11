
#include <cmath>
#include <iostream>
#include "lib.h"
#include <chrono>
using namespace std;
using namespace std::chrono;

double improved_MC(double *);
//     Main function begins here
int main()
{
     int n;
     double x[6], y1, y2, r, fx; //x = [r1,r2,theta1,theta2,phi1,phi2]
     cout << "Read in the number of Monte-Carlo samples"<< endl;
     cin >> n;
     
     time_point<high_resolution_clock> start, end;
     start = high_resolution_clock::now();
     double int_mc = 0.;  double variance = 0.;
     double sum_sigma= 0. ; long idum=-1 ;
     double jacobi_det = 4*pow(acos(-1.),4.)*1/16;

//   evaluate the integral with importance sampling
     for ( int i = 1;  i <= n; i++){
//   x[] contains the random numbers for all dimensions
       //Generate r1 and r2 according to exponential dist.
       y1 = ran0(&idum);
       y2 = ran0(&idum);
       x[0] = -log(1-y1)/4.;   //-log(1-y1)/(2*alpha)
       x[1] = -log(1-y2)/4.;   //-log(1-y2)/(2*alpha)
       for (int j = 2; j< 4; j++) {
           x[j]=acos(-1.)*ran0(&idum);   //pi*random
       }
       for (int j = 4; j< 6; j++) {
           x[j]=2*acos(-1.)*ran0(&idum); //2*pi*random
       }
       fx=improved_MC(x);
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
     cout << "Standard deviation = "<< jacobi_det*sqrt(variance/n) <<  " Integral = " << jacobi_det*int_mc << " exact= " << 5*M_PI*M_PI/(16*16) << " Time = " << time << endl;
     return 0;
}  // end of main program

// this function defines the integrand to integrate

double  improved_MC(double *x)
{
   double alpha = 2.;
// evaluate the different terms of the exponential
   double cosb = cos(x[2])*cos(x[3])+sin(x[2])*sin(x[3])*cos(x[4]-x[5]);
   double deno = sqrt(x[0]*x[0]+x[1]*x[1]-2*x[0]*x[1]*cosb);
   double value = x[0]*x[0]*x[1]*x[1]*sin(x[2])*sin(x[3])/deno;
   return value;

} // end function for the integrand
