
#include <cmath>
#include <iostream>
#include "lib.h"
using namespace std;

double brute_force_MC(double *);
//     Main function begins here     
int main()
{
     int n;
     double x[6], y, fx; 
     printf("Read in the number of Monte-Carlo samples\n");
     scanf("%d", &n);

     double int_mc = 0.;  double variance = 0.;
     double sum_sigma= 0. ; long idum=-1 ;  
     double length=5.; // we fix the max size of the box to L=5
     double volume=pow((2*length),6);

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
//   final output 
     printf("%d standard deviation= %12.5E, Inum= %12.5E, exact= %f", n, volume*sqrt(variance/n), volume*int_mc, 5*M_PI*M_PI/(16*16)); 
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
