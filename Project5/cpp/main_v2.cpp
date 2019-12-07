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
  string filename;
  if (argc < 2)
  {
    cout << "Please read in output file name on the same line." << endl;
    exit(1);
  }
  else
  {
    filename = argv[1];
    filename = "results/" + filename;
  }
  int mcs;
  double omega;
  cout << "Please enter the number of monte carlo cycles" << endl;
  cin >> mcs;
  cout << "Please read in frequency (omega)" << endl;
  cin >> omega;

  ofile.open(filename);
  ofile << " Alpha    E w/o Coloumb   Var w/o Coloumb      E w Coloumb"
           "    V w Coloumb  Accepted moves w/o Coloumb  Accepted moves w Coloumb"
        << endl;

  int accepted_moves_wo_coloumb;
  int accepted_moves_w_coloumb;
  double energy_wo_coloumb, variance_wo_coloumb, r12_wo_coloumb, KE_wo_coloumb, var_KE1;
  double energy_w_coloumb, variance_w_coloumb, r12_w_coloumb, KE_w_coloumb, var_KE2;

  double beta = 0;
  for (double alpha = 0.2; alpha < 1.6; alpha += 0.02)
  {
    var_mc(energy_wo_coloumb, variance_wo_coloumb, r12_wo_coloumb,
           accepted_moves_wo_coloumb, mcs, alpha, beta, omega,
           TrialWaveFunction1, E1, KE_wo_coloumb, var_KE1);
    var_mc(energy_w_coloumb, variance_w_coloumb, r12_w_coloumb,
           accepted_moves_w_coloumb, mcs, alpha, beta, omega,
           TrialWaveFunction1, E_repuls, KE_w_coloumb, var_KE2);

    ofile << setprecision(8) << setw(6) << alpha;
    ofile << setprecision(8) << setw(17) << energy_wo_coloumb;
    ofile << setprecision(8) << setw(18) << variance_wo_coloumb;
    ofile << setprecision(8) << setw(17) << energy_w_coloumb;
    ofile << setprecision(8) << setw(15) << variance_w_coloumb;
    ofile << setprecision(8) << setw(28) << accepted_moves_wo_coloumb;
    ofile << setprecision(8) << setw(26) << accepted_moves_w_coloumb << endl;
  }
  return 0;
}

//Define energy func+
