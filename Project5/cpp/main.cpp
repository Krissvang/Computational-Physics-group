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
  cout << "Please enter the number of monte carlo cylcles" << endl;
  cin >> mcs;

  char *filename = argv[1];
  ofile.open(filename);
  ofile << "Alpha       E w/o Colomb   Var w/o Colomb   E w Colomb   "
        " V w Colomb  Accepted moves" << endl;

  mat r(2, 3);

  double alpha = 1;
  double beta = 0;
  double omega = 1.;

  int accepted_moves = 0;

  for (double alpha = 0.2; alpha < 1.6; alpha += 0.1)
  {
    double h = FindOptimal_h(alpha, beta, omega, mcs/10, TrialWaveFunction1);
    r = init_pos();

    double energy_wo_colomb = 0;
    double energy_w_colomb = 0;

    double energy2_wo_colomb = energy_wo_colomb * energy_wo_colomb;
    double energy2_w_colomb = energy_w_colomb*energy_w_colomb;
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
      
      double w = P1(local_r, alpha, omega) / P1(r, alpha, omega);
      if (unif_dist(gen) <= w)
      {
        r = local_r;
        accepted_moves += 1;
      }
      double local_energy_wo_colomb = E1(r, alpha, omega);
      double local_energy_w_colomb = E1(r, alpha, omega)+E_repuls(r,alpha,omega);

      energy_wo_colomb += local_energy_wo_colomb;
      energy_w_colomb += local_energy_w_colomb;

      energy2_wo_colomb += local_energy_wo_colomb * local_energy_wo_colomb;
      energy2_w_colomb += local_energy_w_colomb * local_energy_w_colomb;

    }
    energy_wo_colomb /= mcs;
    energy_w_colomb /=mcs;

    energy2_wo_colomb /= mcs;
    energy2_w_colomb/=mcs;

    double variance_wo_colomb = abs(energy_wo_colomb * energy_wo_colomb - energy2_wo_colomb);
    double variance_w_colomb = abs(energy_w_colomb * energy_w_colomb - energy2_w_colomb);



    ofile << setprecision(8) << setw(6) << alpha;
    ofile << setprecision(8) << setw(15) << energy_wo_colomb;
    ofile << setprecision(8) << setw(15) << variance_wo_colomb;
    ofile << setprecision(8) << setw(15) << energy_w_colomb;
    ofile << setprecision(8) << setw(15) << variance_w_colomb;
    ofile << setprecision(8) << setw(12) << accepted_moves << endl;
    accepted_moves = 0;
  }
  return 0;
}

//Define energy func+
