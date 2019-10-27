#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "ising.h"
using namespace std;
using namespace std::chrono;

ofstream ofile;

int main(int argc, char const *argv[]) {
  int n_spins, mcs;
  double T_low, T_high, dT;
  string filename, init;
  // read in various parameters from terminal:
  cout << "Read in lattice size (number of spins in each dimension):" << endl;
  cin >> n_spins;
  cout << "Read in number of Monte Carlo cycles:" << endl;
  cin >> mcs;
  cout << "Read in low temperature limit:" << endl;
  cin >> T_low;
  cout << "Read in high temperature limit:" << endl;
  cin >> T_high;
  cout << "Read in step between temperatures:" << endl;
  cin >> dT;
  cout << "Read in spin initialization (r for random, any other "
  "value for ordered):" << endl;
  cin >> init;
  cout << "Read in output file name:" << endl;
  cin >> filename;

  filename = "results/"+filename+".txt";
  ofile.open(filename);
  ofile << n_spins << "x" << n_spins << " spins, " << mcs <<
  " Monte Carlo cycles";
  if(init == "r") ofile << ", random initial state" << endl;
  else ofile << ", ordered initial state" << endl;
  ofile << "Temperature     Avg. energy   Heat capacity  Avg. magnetization   "
  "Avg. magnetization (abs. value)  Susceptibility            Time" << endl;
  time_point<high_resolution_clock> start, end;
  double E_avg, heatcap, M_avg, M_abs_avg, susc, time;
  for(double temperature = T_low; temperature < T_high; temperature += dT){
    start = high_resolution_clock::now();
    ising(n_spins,mcs,temperature,init,E_avg,heatcap,M_avg,M_abs_avg,susc);
    end = high_resolution_clock::now();
    duration<double> elapsed = end-start;
    time = elapsed.count();
    ofile << setprecision(2) << setw(11) << temperature;
    ofile << setprecision(8) << setw(16) << E_avg;
    ofile << setprecision(8) << setw(16) << heatcap;
    ofile << setprecision(8) << setw(20) << M_avg;
    ofile << setprecision(8) << setw(34) << M_abs_avg;
    ofile << setprecision(8) << setw(16) << susc;
    ofile << setprecision(8) << setw(16) << time << endl;
  }
  ofile.close();
  return 0;
}
