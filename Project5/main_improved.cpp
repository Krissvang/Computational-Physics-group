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

int main(int argc, char *argv[])
{
  //Initialize random number generator using mt19937
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> pos_dist(-1.0, 1.0);
  uniform_real_distribution<double> unif_dist(0, 1);
  uniform_real_distribution<double> move_dist(-0.5, 0.5);

  int mcs;
  double omega, beta;
  cout << "Please enter the number of monte carlo cylcles" << endl;
  cin >> mcs;
  cout << "Please enter the frequency (omega)" << endl;
  cin >> omega;
  cout << "Please enter beta" << endl;
  cin >> beta;

  char *filename = argv[1];
  ofile.open(filename);
  ofile << "Alpha       E w/o Colomb   Var w/o Colomb    "
        "  Accepted moves" << endl;

  mat r(2, 3);

  int accepted_moves = 0;

  for (double alpha = 0.2; alpha < 1.6; alpha += 0.1)
  {
    double h = FindOptimal_h(alpha, beta, omega, mcs/10, TrialWaveFunction2);
    r = init_pos();

    double energy = 0;
    double energy2 = energy * energy;
    for (int i = 0; i < mcs; i++)
    {
      mat local_r(2, 3);
      for (int i = 0; i < 2; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          local_r(i,j)=r(i,j)+h*move_dist(gen);
        }
      }

      double w = P2(local_r, alpha, omega, beta) / P2(r, alpha, omega, beta);
      if (unif_dist(gen) <= w)
      {
        r = local_r;
        accepted_moves += 1;
      }
      double local_energy = E2(r, alpha, omega, beta);
      energy += local_energy;
      energy2 += local_energy * local_energy;

    }
    energy /= mcs;
    energy2 /= mcs;

    double variance = abs(energy * energy - energy2);



    ofile << setprecision(8) << setw(6) << alpha;
    ofile << setprecision(8) << setw(15) << energy;
    ofile << setprecision(8) << setw(15) << variance;
    ofile << setprecision(8) << setw(12) << accepted_moves << endl;
    accepted_moves = 0;
  }
  return 0;
}

//Define energy func+
