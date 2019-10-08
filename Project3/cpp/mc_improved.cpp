
#include <cmath>
#include <iostream>
#include "lib.h"
using namespace std;

double improved_MC(double *);
//     Main function begins here     
int main()
{
     int n;
     double x[6], y, fx; 
     printf("Read in the number of Monte-Carlo samples\n");
     scanf("%d", &n);

     double int_mc = 0.;  double variance = 0.;
     double sum_sigma= 0. ; long idum=-1 ;  
     double max_r=1.; // we fix the max size of the box to L=5
     double volume=(4/3)*acos(-1.)*pow(max_r,3);
     double jacobi_det = 4*pow(acos(-1.),4.)*1/16;

//   evaluate the integral with importance sampling    
     for ( int i = 1;  i <= n; i++){
//   x[] contains the random numbers for all dimensions
       for (int j = 0; j< 2; j++) {
           y=ran0(&idum)*(1-exp(-max_r));
           x[j]=-log(1.-y);
           cout << x[j] << endl;
           
       }
       for (int j = 2; j< 4; j++) {
           x[j]=2*acos(-1.)*ran0(&idum);  //2*pi*random
       }
       for (int j = 4; j< 6; j++) {
           x[j]=acos(-1.)*ran0(&idum);    
       }
       fx=improved_MC(x); 
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
 
double  improved_MC(double *x) 
{
   double alpha = 2.;
// evaluate the different terms of the exponential
   double exp1 = - alpha * (x[0]+x[1]) ;
   double cosb = cos(x[2])*cos(x[3])+sin(x[2])*sin(x[3])*cos(x[4]-x[5]);
   double deno = sqrt(x[0]*x[0]+x[1]*x[1]-2*x[0]*x[1]*cosb);
   double value = x[0]*x[0]*x[1]*x[1]*sin(x[2])*sin(x[3])*exp(exp1)/deno;
   return value;
   
} // end function for the integrand
