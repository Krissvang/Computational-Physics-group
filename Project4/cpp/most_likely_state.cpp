#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "ising.h"
using namespace std;
using namespace std::chrono;

ofstream ofile;

int main(int argc, char const *argv[])
{
  int n_spins, mcs, steady_start;
  double mcs_start, mcs_end, dM, temperature;
  string filename, init;
  // read in various parameters from terminal:
  cout << "Read in lattice size (number of spins in each dimension):" << endl;
  cin >> n_spins;
  cout << "Read in when to start Monete carlo cycles:" << endl;
  cin >> steady_start;
  cout << "Read in number of Monte Carlo cycles start:" << endl;
  cin >> mcs_start;
  cout << "Read in number of Monte Carlo cycles end:" << endl;
  cin >> mcs_end;
  cout << "Read in step between Monte Carlo cycles:" << endl;
  cin >> dM;
  cout << "Read in temperature:" << endl;
  cin >> temperature;
  cout << "Read in spin initialization (r for random, any other "
          "value for ordered):"
       << endl;
  cin >> init;
  cout << "Read in output file name:" << endl;
  cin >> filename;

  filename = "results/" + filename + ".txt";
  ofile.open(filename);
  ofile << n_spins << "x" << n_spins << " spins, " << temperature << " Temperature";
  if (init == "r")
    ofile << ", random initial state" << endl;
  else
    ofile << ", ordered initial state" << endl;
  ofile << "Monte Carlo cycles     Avg. energy   Heat capacity  Avg. magnetization   "
           "Avg. magnetization (abs. value)  Susceptibility            Time    Accepted configs"
        << endl;
  time_point<high_resolution_clock> start, end;
  double E_avg, heatcap, M_avg, M_abs_avg, susc, time;
  vector<int> count(n_spins * n_spins + 1);
  for (double mcs = mcs_start; mcs < mcs_end; mcs += dM)
  {
    int count_configs = 0;
    start = high_resolution_clock::now();
    ising(n_spins, mcs, temperature, init, E_avg, heatcap, M_avg, M_abs_avg, susc, count, steady_start, count_configs);
    end = high_resolution_clock::now();
    duration<double> elapsed = end - start;
    time = elapsed.count();
    ofile << setprecision(8) << setw(19) << mcs;
    ofile << setprecision(8) << setw(16) << E_avg;
    ofile << setprecision(8) << setw(16) << heatcap;
    ofile << setprecision(8) << setw(20) << M_avg;
    ofile << setprecision(8) << setw(34) << M_abs_avg;
    ofile << setprecision(8) << setw(16) << susc;
    ofile << setprecision(8) << setw(16) << time;
    ofile << setprecision(8) << setw(16) << count_configs << endl;
  }
  ofile.close();
  return 0;
}
