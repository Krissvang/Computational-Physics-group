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
  double omega, beta;
  cout << "Please enter the number of monte carlo cylcles" << endl;
  cin >> mcs;
  cout << "Please enter the frequency (omega)" << endl;
  cin >> omega;
  cout << "Please enter beta" << endl;
  cin >> beta;

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
  ofile << " Alpha              E            Var            r12"
           "  Accepted moves"
        << endl;

  int accepted_moves;
  double energy, variance, r12;
  for (double alpha = 0.2; alpha < 1.6; alpha += 0.1)
  {
    var_mc(energy, variance, r12, accepted_moves, mcs, alpha, beta, omega, TrialWaveFunction1, E_repuls);

    ofile << setprecision(8) << setw(6) << alpha;
    ofile << setprecision(8) << setw(15) << energy;
    ofile << setprecision(8) << setw(15) << variance;
    ofile << setprecision(8) << setw(15) << r12;
    ofile << setprecision(8) << setw(16) << accepted_moves << endl;
  }
  return 0;
}
