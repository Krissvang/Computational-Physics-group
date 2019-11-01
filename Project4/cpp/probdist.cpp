#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "ising.h"
#include <mpi.h>
#include <vector>
using namespace std;
using namespace std::chrono;

ofstream ofile;

int main(int argc, char *argv[]) {
  int numprocs, my_rank;
  //Starts the parallizing
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
  int n_spins, mcs, steady_start;
  double T;
  string filename, init;
  if(argc < 7){
    if(my_rank == 0){
      cout << "Please read in the following on the same line, in this order: "
      "Number of spins in each direction, #Monte Carlo cycles, temperature, "
      "#Monte Carlo cycles before equilibrium is reached, spin initialization "
      "(r for random, any other value for ordered), and output file name." <<
      endl;
    }
    exit(1);
  }
  else{
    n_spins = atoi(argv[1]); mcs = atoi(argv[2]); T = atof(argv[3]);
    steady_start = atoi(argv[4]); init = argv[5]; filename = argv[6];
  }

  filename = "results/"+filename+".txt";
  if(my_rank == 1){
    ofile.open(filename);
    ofile << n_spins << "x" << n_spins << " spins, " << mcs <<
    " Monte Carlo cycles, temperature = " << T << ", " << numprocs <<
    " processes, started collecting data after " << steady_start << " cycles.";
    if(init == "r") ofile << ", random initial state" << endl;
    else ofile << ", ordered initial state" << endl;
    ofile << "Energy per spin:     Count:" << endl;
  }
  double E_avg, heatcap, M_avg, M_abs_avg, susc;
  double local_E_avg, local_heatcap, local_M_avg, local_M_abs_avg, local_susc;
  int local_mcs = mcs/numprocs;
  vector<int> local_count(n_spins*n_spins+1);
  int count_i, local_count_i;
  double E_min = -2*n_spins*n_spins;
  ising(n_spins,local_mcs,T,init,local_E_avg,local_heatcap,
    local_M_avg,local_M_abs_avg,local_susc, local_count, steady_start);
  for(int i=0; i < local_count.size(); i++){
    local_count_i = local_count[i];
    MPI_Reduce(&local_count_i, &count_i, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(my_rank == 0){
      double E_i =  ((float) (E_min+4*i))/(n_spins*n_spins);
      ofile << setprecision(8) << setw(16) << E_i << setw(11) << count_i
      << endl;
    }
  }

  MPI_Finalize();
  ofile.close();
  return 0;
}
