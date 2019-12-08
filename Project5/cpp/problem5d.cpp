#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <string>
#include <vector>
#include "varmontecarlo.h"

using namespace arma;
using namespace std;

ofstream ofile;

int main(int argc, char *argv[])
{
  int mcs;
  double omega, fixed_par;
  string fixed_par_name;
  cout << "Please enter the number of monte carlo cylcles" << endl;
  cin >> mcs;
  cout << "Please enter the frequency (omega)" << endl;
  cin >> omega;
  cout << "Please enter the parameter to fix (alpha or beta)" << endl;
  cin >> fixed_par_name;
  cout << "Enter its value" << endl;
  cin >> fixed_par;

  string filename;
  if (argc < 2)
  {
    cout << "Please read in filename on the same line." << endl;
    exit(1);
  }
  else
  {
    filename = argv[1];
    filename = "results/" + filename;
  }
  ofile.open(filename);

  int accepted_moves;
  double energy, variance, r12, KE, var_KE;

  double parameter_min = 0.1;
  double parameter_max = 1;
  double parameter_step = 0.01;

  if (fixed_par_name == "alpha")
  {
    ofile << "Alpha = " << fixed_par << endl;
    ofile << " Beta              PE            PE Var            KE"
             "            KE Var            r12"
             "     Accepted moves"
          << endl;
    for (double beta = parameter_min; beta < parameter_max; beta += parameter_step)
    {
      var_mc(energy, variance, r12, accepted_moves, mcs, fixed_par, beta, omega, TrialWaveFunction2, E2, KE, var_KE);

      ofile << setprecision(8) << setw(6) << beta;
      ofile << setprecision(8) << setw(15) << energy;
      ofile << setprecision(8) << setw(15) << abs(variance);
      ofile << setprecision(8) << setw(15) << KE;
      ofile << setprecision(8) << setw(15) << var_KE;
      ofile << setprecision(8) << setw(15) << r12;
      ofile << setprecision(8) << setw(16) << accepted_moves << endl;
    }
  }
  if (fixed_par_name == "beta")
  {
    ofile << "Beta = " << fixed_par << endl;
    ofile << " Alpha             PE        PE Var            KE"
             "          KE Var        r12"
             "      Accepted moves"
          << endl;
    for (double alpha = parameter_min; alpha < parameter_max; alpha += parameter_step)
    {
      var_mc(energy, variance, r12, accepted_moves, mcs, alpha, fixed_par, omega, TrialWaveFunction2, E2, KE, var_KE);

      ofile << setprecision(8) << setw(6) << alpha;
      ofile << setprecision(8) << setw(15) << energy;
      ofile << setprecision(8) << setw(15) << abs(variance);
      ofile << setprecision(8) << setw(15) << KE;
      ofile << setprecision(8) << setw(15) << var_KE;
      ofile << setprecision(8) << setw(15) << r12;
      ofile << setprecision(8) << setw(16) << accepted_moves << endl;
    }
  }

  return 0;
}
