
#include <cmath>
#include <iostream>
#include "lib.h"
#include <chrono>
#include <mpi.h>
#include <time.h>
using namespace std;
using namespace std::chrono;

double improved_MC(double *);
//     Main function begins here
int main(int nargs, char* args[])
{
     time_point<high_resolution_clock> start, end;
     time_point<system_clock> time2;
     start = high_resolution_clock::now();
     time2 = system_clock::now();
     duration<double> duration_in_seconds =duration<double>(time2.time_since_epoch());
     long t2= duration_in_seconds.count();
     int numprocs, my_rank, i, n;
     n=10000000;
     MPI_Init (&nargs, &args);
     MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
     MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
     int local_n=n/numprocs;
     double x[6], y1, y2, r, fx; //x = [r1,r2,theta1,theta2,phi1,phi2]

     double int_mc = 0.;  double variance = 0.;
     double local_int_mc=0.; double local_variance=0.;
     double sum_sigma= 0. ; long idum=t2+my_rank;
     cout << t2*my_rank << endl;
     double local_sum_sigma=0.;
     double jacobi_det = 4*pow(acos(-1.),4.)*1/16;

//   evaluate the integral with importance sampling
     for ( int i = 1;  i <= local_n; i++){
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
       local_int_mc += fx;
       local_sum_sigma += fx*fx;
     }
     local_int_mc = local_int_mc;
     local_sum_sigma = local_sum_sigma;
     local_variance=local_sum_sigma-local_int_mc*local_int_mc;
     MPI_Reduce(&local_int_mc, &int_mc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     MPI_Reduce(&local_sum_sigma, &sum_sigma, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     int_mc*=1./n;
     sum_sigma*=1./n;
     variance = sum_sigma-int_mc*int_mc;
     end = high_resolution_clock::now();
     duration<double> elapsed = end-start;
     double time = elapsed.count();
    //   final output
     if( my_rank==0){
     cout << "Standard deviation = "<< jacobi_det*sqrt(variance/n) <<  " Integral = " << jacobi_det*int_mc << " exact= " << 5*M_PI*M_PI/(16*16) << " Time = " << time << endl; 
     }
     MPI_Finalize();
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
