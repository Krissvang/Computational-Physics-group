#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <iomanip>
#include <armadillo>
#include <vector>
#include "lib.h"

using namespace arma;
using namespace std;


ofstream ofile;

int main(int argc, char *argv[]){
  int mcs;
  cout << "Please enter the number of monte carlo cylcles" << endl;
  cin >> mcs;

  char *filename = argv[1];
  ofile.open(filename);
  ofile << "Alpha       Beta       Omega     Stepsize" << endl;

  mat r(2, 3);
  r = init_pos();

  double beta = 0;
  double omega = 1.;


  for (double alpha = 0.2; alpha < 1.6; alpha += 0.1)
  {
    double h = FindOptimal_h(alpha, beta, omega, mcs, TrialWaveFunction1);
    ofile << setprecision(8) << setw(4) << alpha;
    ofile << setprecision(8) << setw(12) << beta;
    ofile << setprecision(8) << setw(12) << omega;
    ofile << setprecision(8) << setw(12) << h << endl;
  }

  return 0;
}