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



int main{

  
  int mcs;
  cout << "Please enter the number of monte carlo cylcles" << endl;
  cin >> mcs;

  char *filename = argv[1];
  ofile.open(filename);
  ofile << "Alpha       Beta       Omega     Stepsize" << endl;

  mat r(2, 3);
  r = init_pos();

  double alpha = 1;
  double beta = 0;
  double omega = 1.;

  int accepted_moves = 0;

  for (double alpha = 0.2; alpha < 1.6; alpha += 0.1)
  {
    double h = FindOptimal_h(alpha, beta, omega, mcs, TrialWaveFunction1);
  }
    double variance = abs(energy * energy - energy2);
    ofile << setprecision(8) << setw(12) << alpha;
    ofile << setprecision(8) << setw(12) << energy;
    ofile << setprecision(8) << setw(12) << variance;
    ofile << setprecision(8) << setw(12) << accepted_moves << endl;
    accepted_moves = 0;
  }




  return 0;
}