
#include <cmath>
#include <iostream>
#include <iomanip>
#include "lib.h"
#include <chrono>
#include <mpi.h>
#include <time.h>
#include "montecarlo.h"
using namespace std;
using namespace std::chrono;

double improved_MC(double *);
//     Main function begins here
int main(int nargs, char* args[])
{   
    //Takes the time
    time_point<high_resolution_clock> start, end;
    start = high_resolution_clock::now();

    int n=pow(10,atoi(args[1]));
    int numprocs, my_rank, i;
    //Starts the parallizing
    MPI_Init (&nargs, &args);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    time_point<system_clock> time2;
    time2 = system_clock::now();
    duration<double> duration_in_seconds =duration<double>(time2.time_since_epoch());
    long t2= duration_in_seconds.count();
    
    //Needed loocal and global variables
    int local_n=n/numprocs;
    double int_mc = 0.;  double variance = 0.;
    double local_int_mc=0.; double local_variance=0.;
    double sum_sigma= 0. ; long idum=t2+my_rank;
    double local_sum_sigma=0.;
    double jacobi_det = 4*pow(acos(-1.),4.)*1/16;
    double time=0.; double local_std_dev;
    //runs the algorithm
    mc_improved(&improved_MC,local_n,local_int_mc,local_std_dev,time,local_sum_sigma,t2);
    //Fixes the dimentions 
    local_int_mc*=local_n;
    local_sum_sigma*=local_n;
    //Sums the computed mean values together
    MPI_Reduce(&local_int_mc, &int_mc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sum_sigma, &sum_sigma, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //Fixes the dimentions 
    int_mc*=1./n/jacobi_det;
    sum_sigma*=1./n;
    variance = sum_sigma-int_mc*int_mc;
    
    //   final output
    if( my_rank==0){
    end = high_resolution_clock::now();
    duration<double> elapsed = end-start;
    time = elapsed.count();
    cout << "Standard deviation = "<< jacobi_det*sqrt(variance/n) <<  " Integral = " << setprecision(10) << jacobi_det*int_mc << " exact= " << setprecision(10) << 5*M_PI*M_PI/(16*16) << " Time = " << time << endl; 
    }
    MPI_Finalize();
    return 0;
}  // end of main program

// this function defines the integrand to integrate
double  improved_MC(double *x)
{
   double cosb = cos(x[2])*cos(x[3])+sin(x[2])*sin(x[3])*cos(x[4]-x[5]);
   double deno = sqrt(x[0]*x[0]+x[1]*x[1]-2*x[0]*x[1]*cosb);
   double value = x[0]*x[0]*x[1]*x[1]*sin(x[2])*sin(x[3])/deno;
   return value;

} // end function for the integrand
