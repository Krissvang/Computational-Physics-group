#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include "varmontecarlo.h"

using namespace arma;
using namespace std;

ofstream ofile;

int main(int argc, char *argv[])
{
  int mcs;
  double beta;
  cout << "Please enter the number of monte carlo cylcles" << endl;
  cin >> mcs;

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

  ofile << " Omega          E        Var E              KE"
           "          KE Var           PE wo C     Var PE wo C   "
           "   PE w C        Var PE w C         r12"
           "          Accepted moves"
        << endl;

  int accepted_moves;
  double energy, variance, r12, KE, var_KE, PE_wo, PE_w, var_PE_w, var_PE_wo;
  double alpha = 0.87;
  ofile << "Alpha = " << alpha << endl;
  for (double omega = 0.2; omega < 1.5; omega += 0.05)
  {
    var_mc(energy, variance, r12, accepted_moves, mcs, alpha, beta, omega, TrialWaveFunction1, E_repuls, KE, var_KE, PE_wo, var_PE_wo, PE_w, var_PE_w);

    ofile << setprecision(8) << setw(6) << omega;
    ofile << setprecision(8) << setw(15) << energy;
    ofile << setprecision(8) << setw(15) << abs(variance);
    ofile << setprecision(8) << setw(15) << KE;
    ofile << setprecision(8) << setw(15) << var_KE;
    ofile << setprecision(8) << setw(15) << PE_wo;
    ofile << setprecision(8) << setw(15) << var_PE_wo;
    ofile << setprecision(8) << setw(15) << PE_w;
    ofile << setprecision(8) << setw(15) << var_PE_w;
    ofile << setprecision(8) << setw(15) << r12;
    ofile << setprecision(8) << setw(16) << accepted_moves << endl;
  }
  return 0;
}
