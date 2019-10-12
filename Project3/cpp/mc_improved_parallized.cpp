
#include <cmath>
#include <iostream>
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
     time_point<high_resolution_clock> start, end;
     time_point<system_clock> time2;
     start = high_resolution_clock::now();
     time2 = system_clock::now();
     duration<double> duration_in_seconds =duration<double>(time2.time_since_epoch());
     long t2= duration_in_seconds.count();
     int numprocs, my_rank, i;
     MPI_Init (&nargs, &args);
     MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
     MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
     int n;
     if(nargs < 2){
       if(my_rank == 0){
         cout << "Please read in number of Monte Carlo "
         "samples on the same line." << endl;
       }
       exit(1);
     }
     else n = atoi(args[1]);
     int local_n=n/numprocs;
     double int_mc, std_dev, local_int_mc, local_std_dev, local_time;
     mc_improved(&improved_MC,local_n,local_int_mc,local_std_dev,local_time);
     MPI_Reduce(&local_int_mc, &int_mc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     MPI_Reduce(&local_std_dev, &std_dev, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     int_mc/=numprocs;
     std_dev/=numprocs*sqrt(numprocs);
     end = high_resolution_clock::now();
     duration<double> elapsed = end-start;
     double time = elapsed.count();
    //   final output
     if( my_rank==0){
     cout << "Standard deviation = "<< std_dev <<  " Integral = " << int_mc <<
     " exact= " << 5*M_PI*M_PI/(16*16) << " Time = " << time << endl;
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
