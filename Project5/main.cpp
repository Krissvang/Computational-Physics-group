#include <iostream>
#include <cmath>
#include <random>

using namespace std;

double E1(double r1, double r2, double alph, double w);
double P1(double r1, double r2, double alph, double w);

int main()
{
  //Initialize random number generator using mt19937
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> RNG(0.0, 1.0);
  int mcs;
  cout << "Please enter the number of monte carlo cylcles" << endl;
  cin >> mcs;

  double r1 = 0;
  double r2 = 0;
  double alph = 1.;
  double omega = 1.;
  double delta = .07;
  double cutoff = 10;
  int accepted_moves = 0;

  double energy = E1(r1, r2, alph, omega);
  double energy2 = energy * energy;
  for (int i = 0; i < mcs; i++)
  {
    double x1 = RNG(gen);
    double x2 = RNG(gen);
    x1 *= cutoff;
    x2 *= cutoff;
    double local_r1 = r1 + delta * x1;
    double local_r2 = r2 + delta * x2;

    double local_energy = E1(local_r1, local_r2, alph, omega);
    double w = P1(local_r1, local_r2, alph, omega) / P1(r1, r2, alph, omega);
    if (RNG(gen) <= w)
    {
      r1 = local_r1;
      r2 = local_r2;
      accepted_moves += 1;
    }
    energy += local_energy;
    energy2 += local_energy * local_energy;
  }
  cout << "Energy: " << energy / mcs << "   Variance: " << energy2 / mcs << endl;
  return 0;
}

double E1(double r1, double r2, double alph, double w)
{
  double E = 0.5 * w * w * (r1 * r1 + r2 * r2) * (1 - alph * alph) + 3 * alph * w;
  return E;
}

double P1(double r1, double r2, double alph, double w)
{
  double P = abs(alph * w / acos(-1)) * exp(-alph * w * (r1 * r1 + r2 * r2));
  return P;
}

//Define energy func+
