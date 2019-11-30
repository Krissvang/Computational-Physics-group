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

int main(int argc, char *argv[])
{
  mat r(2, 3);
  r.zeros();
  r = init_pos();
  int mcs = 1000;
  double alpha = 1;
  double omega = 1;
  double beta = 0;
  double h = FindOptimal_h(alpha, beta, omega, mcs, TrialWaveFunction1);
  double r1 = r_squared(r);
  cout << r << endl;
  return 0;
}
