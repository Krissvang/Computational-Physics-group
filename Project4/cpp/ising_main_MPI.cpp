#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "ising.h"
#include <mpi.h>
using namespace std;
using namespace std::chrono;

ofstream ofile;

int main(int argc, char *argv[]) {
  int numprocs, my_rank;
  int n_spins, mcs, steady_start;
  double T_low, T_high, dT;
  string filename, init;
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  if(argc < 8){
   if(my_rank == 0){
    cout << "Please read in the following on the same line, in this order: "
    "Number of spins in each direction, #Monte Carlo cycles, "
    "#Monte Carlo cycles before equilibrium is reached, low temperature limit temperature, high temperature limit, step between temperatures , spin initialization "
    "(r for random, any other value for ordered), and output file name." <<
    endl;
    }
  exit(1);
  }
  // read in various parameters from terminal:
  else{
  n_spins=atoi(argv[1]); mcs=atoi(argv[2]); steady_start=atoi(argv[3]); T_low=atoi(argv[4]); T_high=atoi(argv[5]); dT=atoi(argv[6]); init=atoi(argv[7]); filename=atoi(argv[8]);
  }

  filename = "results/"+filename+".txt";
  if(my_rank==0){
    ofile.open(filename);
    ofile << n_spins << "x" << n_spins << " spins, " << mcs <<
    " Monte Carlo cycles";
    if(init == "r") ofile << ", random initial state" << endl;
    else ofile << ", ordered initial state" << endl;
    ofile << "Temperature     Avg. energy   Heat capacity  Avg. magnetization   "
    "Avg. magnetization (abs. value)  Susceptibility            Time    " << endl;
  }
  time_point<high_resolution_clock> start, end;
  double E_avg, heatcap, M_avg, M_abs_avg, susc, time;
  double local_E_avg, local_heatcap, local_M_avg, local_M_abs_avg, local_susc;
  int local_mcs=mcs/numprocs;
  int local_steady_start = steady_start/numprocs;
  int count_configs =0;
  vector<int> local_count(n_spins*n_spins+1);
  for(double temperature = T_low; temperature < T_high; temperature += dT){

    start = high_resolution_clock::now();
    ising(n_spins,local_mcs,temperature,init,local_E_avg,local_heatcap,local_M_avg,local_M_abs_avg,local_susc,local_count,local_steady_start,count_configs);
    end = high_resolution_clock::now();
    duration<double> elapsed = end-start;
    time = elapsed.count();

    MPI_Reduce(&local_E_avg, &E_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&local_heatcap, &heatcap, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&local_M_avg, &M_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&local_M_abs_avg, &M_abs_avg, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&local_susc, &susc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
   



    if(my_rank==0){
      E_avg/=numprocs;
      heatcap/=numprocs;
      M_avg/=numprocs;
      M_abs_avg/=numprocs;
      susc/=numprocs;
      cout << temperature << endl;
      ofile << setprecision(2) << setw(11) << temperature;
      ofile << setprecision(8) << setw(16) << E_avg;
      ofile << setprecision(8) << setw(16) << heatcap;
      ofile << setprecision(8) << setw(20) << M_avg;
      ofile << setprecision(8) << setw(34) << M_abs_avg;
      ofile << setprecision(8) << setw(16) << susc;
      ofile << setprecision(8) << setw(16) << time << endl;
    }
  }
  MPI_Finalize();
  ofile.close();
  return 0;
}
