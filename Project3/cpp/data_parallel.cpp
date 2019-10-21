#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include "lib.h"
#include "montecarlo.h"
#include <string>
#include <chrono>

using namespace std::chrono;
using namespace std;

double  improved_MC(double *x);

ofstream ofile;

int main(int argc, char *argv[]) {
  //Takes the time
  time_point<high_resolution_clock> start, end;

  int numprocs, my_rank, i;
  //Starts the parallizing
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  int n;
  string filename;
  if(argc < 3){
    if(my_rank == 0){
      cout << "Please read in number of samples and filename on the same line" << endl;
    }
    exit(1);
  }
  else{
    n = atoi(argv[1]);
    filename = argv[2];
  }
  filename = "mc_results/parallel_data_"+filename+"_"+to_string(numprocs)+"_"+to_string(n)+".txt";
  if(my_rank == 0) ofile.open(filename);

  n=pow(10,n);

  time_point<system_clock> time2;
  time2 = system_clock::now();
  duration<double> duration_in_seconds =duration<double>(time2.time_since_epoch());
  long t2= duration_in_seconds.count();
  t2 += my_rank;

  //Needed loocal and global variables
  int local_n=n/numprocs;
  double int_mc = 0.;  double variance = 0.;
  double local_int_mc=0.;
  double sum_f2= 0.;
  double local_sum_f2=0.;
  double jacobi_det = 4*pow(acos(-1.),4.)*1/16;
  double time=0.; double local_std_dev;

  for (int i = 0; i < 10; i++)
  {
    start = high_resolution_clock::now();
    //runs the algorithm
    mc_improved(&improved_MC,local_n,local_int_mc,local_std_dev,time,local_sum_f2,t2+42*i);
    //Recover local sums of f and f^2
    local_int_mc*=local_n;
    local_sum_f2*=local_n;
    //Sums the computed local sums together
    MPI_Reduce(&local_int_mc, &int_mc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_sum_f2, &sum_f2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //Find global mean of f and f^2
    int_mc*=1./n/jacobi_det;    //=<f> for all processes
    sum_f2*=1./n;               //=<f^2> for all processes
    variance = sum_f2-int_mc*int_mc;

    // final output
    if( my_rank==0){
    end = high_resolution_clock::now();
    duration<double> elapsed = end-start;
    time = elapsed.count();
    ofile << setprecision(10)<< jacobi_det*sqrt(variance/n) <<  setw(30) << setprecision(20) << jacobi_det*int_mc << setprecision(10) << setw(20)<< time << endl;
    }
  }

  MPI_Finalize();
  ofile.close();
  return 0;
}

double  improved_MC(double *x)
{
   double cosb = cos(x[2])*cos(x[3])+sin(x[2])*sin(x[3])*cos(x[4]-x[5]);
   double deno = sqrt(x[0]*x[0]+x[1]*x[1]-2*x[0]*x[1]*cosb);
   double value = x[0]*x[0]*x[1]*x[1]*sin(x[2])*sin(x[3])/deno;
   return value;

} // end function for the integrand
