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

// This program lets the user pick which parameter
// to fix, and let the other one run.
// It then writes the energies, variance
// and all other needed parameters to a file.

int main(int argc, char *argv[])
{
  //Needed parameters
  int mcs;
  double omega, fixed_par;
  string fixed_par_name;

  //Inputs
  cout << "Please enter the number of monte carlo cylcles" << endl;
  cin >> mcs;
  cout << "Please enter the frequency (omega)" << endl;
  cin >> omega;
  cout << "Please enter the parameter to fix (alpha, beta or both)" << endl;
  cin >> fixed_par_name;
  mcs = pow(10, mcs);
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
  //Needed variables
  int accepted_moves;
  double energy, variance, r12, KE, var_KE, PE_wo_C, var_PE_wo,
      PE_w_C, var_PE_w;

  //Decides the interval and stepsize to run.
  //for Alpha
  double parameter_min_a = 0.99;
  double parameter_max_a = 1.005;
  double parameter_step_a = 0.001;
  //For Beta
  double parameter_min_b = 0.27;
  double parameter_max_b = 0.295;
  double parameter_step_b = 0.001;

  //Runs for fixed alpha
  if (fixed_par_name == "alpha")
  {
    cout << "Enter its value" << endl;
    cin >> fixed_par;
    ofile << "Alpha = " << fixed_par << endl;
    ofile << " Beta              E            E Var            KE"
             "            KE Var            r12"
             "     Accepted moves"
          << endl;
    for (double beta = parameter_min_b;
         beta < parameter_max_b; beta += parameter_step_b)
    {
      //runs VMC
      var_mc(energy, variance, r12, accepted_moves, mcs, fixed_par,
             beta, omega, TrialWaveFunction2, E2,
             KE, var_KE, PE_wo_C, var_PE_wo, PE_w_C, var_PE_w);

      //Writes to file
      ofile << setprecision(8) << setw(6) << beta;
      ofile << setprecision(8) << setw(15) << energy;
      ofile << setprecision(8) << setw(15) << abs(variance);
      ofile << setprecision(8) << setw(15) << KE;
      ofile << setprecision(8) << setw(15) << var_KE;
      ofile << setprecision(8) << setw(15) << r12;
      ofile << setprecision(8) << setw(16) << accepted_moves << endl;
    }
  }
  //Runs for fixed beta
  else if (fixed_par_name == "beta")
  {
    cout << "Enter its value" << endl;
    cin >> fixed_par;
    ofile << "Beta = " << fixed_par << endl;
    ofile << " Alpha             E        E Var            KE"
             "          KE Var        r12"
             "      Accepted moves"
          << endl;
    for (double alpha = parameter_min_a;
         alpha < parameter_max_a; alpha += parameter_step_a)
    {
      var_mc(energy, variance, r12, accepted_moves, mcs, alpha,
             fixed_par, omega, TrialWaveFunction2, E2, KE,
             var_KE, PE_wo_C, var_PE_wo, PE_w_C, var_PE_w);

      ofile << setprecision(8) << setw(6) << alpha;
      ofile << setprecision(8) << setw(15) << energy;
      ofile << setprecision(8) << setw(15) << abs(variance);
      ofile << setprecision(8) << setw(15) << KE;
      ofile << setprecision(8) << setw(15) << var_KE;
      ofile << setprecision(8) << setw(15) << r12;
      ofile << setprecision(8) << setw(16) << accepted_moves << endl;
    }
  }
  //Fixes alpha and beta to the found parameters and let omega run
  else if (fixed_par_name == "both")
  {
    double alpha = 1.002;
    double beta = 0.276;
    ofile << "Alpha = " << alpha << ". And Beta = " << beta << endl;
    ofile << " Omega          E        Var E              KE"
             "          KE Var           PE wo C     Var PE wo C   "
             "   PE w C        Var PE w C         r12"
             "          Accepted moves"
          << endl;
    for (double omega = 0.2; omega < 1.05; omega += 0.05)
    {
      var_mc(energy, variance, r12, accepted_moves, mcs, alpha,
             beta, omega, TrialWaveFunction2, E2, KE, var_KE,
             PE_wo_C, var_PE_wo, PE_w_C, var_PE_w);

      ofile << setprecision(8) << setw(6) << omega;
      ofile << setprecision(8) << setw(15) << energy;
      ofile << setprecision(8) << setw(15) << abs(variance);
      ofile << setprecision(8) << setw(15) << KE;
      ofile << setprecision(8) << setw(15) << var_KE;
      ofile << setprecision(8) << setw(15) << PE_wo_C;
      ofile << setprecision(8) << setw(15) << var_PE_wo;
      ofile << setprecision(8) << setw(15) << PE_w_C;
      ofile << setprecision(8) << setw(15) << var_PE_w;
      ofile << setprecision(8) << setw(15) << r12;
      ofile << setprecision(8) << setw(16) << accepted_moves << endl;
    }
  }
  else
  {
    cout << "Please pick either alpha,"
            " beta or both as parameters to fix"
         << endl;
  }

  return 0;
}