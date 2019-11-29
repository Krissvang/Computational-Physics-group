#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <armadillo>
#include <iomanip>
#include "lib.h"

using namespace std;
using namespace arma;

double r_12(mat &r)
{
  double r_12 = 0;
  for (int i = 0; i < 3; i++)
  {
    r_12 += (r(0, i) - r(1, i)) * (r(0, i) - r(1, i));
  }
  return sqrt(r_12);
}

double radi(int i, mat &r)
{
  double ri = 0;
  for (int j = 0; j < 3; j++)
  {
    ri += r(i, j) * r(i, j);
  }
  return sqrt(ri);
}

double r_squared(mat &r)
{
  double r2 = 0;
  for (int i = 0; i < 3; i++)
  {
    r2 += r(0, i) * r(0, i) + r(1, i) * r(1, i);
  }
  return r2;
}

double TrialWaveFunction1(mat &r, double alpha, double beta, double omega)
{
  double argument, wavefunction, r_i;
  argument = wavefunction = 0;

  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      r_i += r(i, j) * r(i, j);
    }
    argument += r_i;
  }
  wavefunction = exp(-0.5 * alpha * omega * argument);
  return wavefunction;
}

double TrialWaveFunction2(mat &r, double alpha, double beta, double omega)
{
  double arg = 1.0 + beta * r_12(r);
  double part2 = exp(r_12(r) / (2 * arg));
  double part1 = TrialWaveFunction1(r, alpha, beta, omega);
  return part1 * part2;
}

mat init_pos()
{
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> pos_dist(-1.0, 1.0);
  mat r(2, 3);
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      r(i, j) = pos_dist(gen);
    }
  }
  return r;
}

//Program that finds optimal step length h
double FindOptimal_h(double alpha, double beta, double omega, int mcs,
                     double (*TrialWaveFunction)(mat &, double, double, double))
{
  random_device rd;
  mt19937_64 gen(rd());
  uniform_real_distribution<double> pos_dist(-1.0, 1.0);
  uniform_real_distribution<double> unif_dist(0, 1);
  uniform_real_distribution<double> move_dist(-0.5, 0.5);

  double h = 10;
  double error = 0.05; //Five percent
  int accepted = 0;    //To store number of accepted transitions
  double acceptance_ratio = accepted / ((double)mcs);
  double w = 0;
  double new_wave_func, old_wave_func;
  bool condition = (0.5 - error) < acceptance_ratio && acceptance_ratio < (0.5 + error);

  while (condition != 1 && h > 0)
  {
    mat r_new(2, 3);
    r_new.zeros();
    mat r_old(2, 3);

    for (int cycle = 0; cycle < mcs; cycle++)
    {

      r_old = init_pos();

      old_wave_func = TrialWaveFunction(r_old, alpha, beta, omega);
      for (int i = 0; i < 2; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          r_new(i, j) = r_old(i, j) + h * move_dist(gen);
        }
      }

      new_wave_func = TrialWaveFunction(r_new, alpha, beta, omega);
      //Metropolis

      w = (new_wave_func * new_wave_func) / (old_wave_func * old_wave_func);
      if (w >= unif_dist(gen))
      {
        r_old = r_new;
        accepted += 1;
        old_wave_func = new_wave_func;
      }
    }

    acceptance_ratio = accepted / ((double)mcs);
    accepted = 0;
    h -= 0.01;
    condition = acceptance_ratio > 0.5 + error;
  }
  return h;
}

double E1(mat &r, double alpha, double omega, double (*r_squared)(mat &))
{
  double E = 0.5 * omega * omega * (r_squared(r)) * (1 - alpha * alpha) + 3 * alpha * omega;
  return E;
}

double E_repuls(mat &r, double alpha, double omega, double (*r_squared)(mat &))
{
  double E2 = 1 / (r_12(r));
  return E2;
}

double P1(mat &r, double alpha, double omega, double (*r_squared)(mat &))
{
  double P = abs(alpha * omega / acos(-1)) * exp(-alpha * omega * (r_squared(r)));
  return P;
}

double Kinetic_E(mat &r, double f, double h, double alpha, double beta, double omega, int mcs,
                 double (*trialFunction)(mat &r, double, double, double))
{
  double f_minus, f_plus, T_Local;
  double hh = h * h;
  mat r_plus(2, 3), r_minus(2, 3);
  r_plus = r_minus = r;
  double h_der = 1 / mcs;    //Step-size for the differentiation
  double h2_der = mcs * mcs; //Inverse squared Step-size for the differentiation

  T_Local = 0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      r_plus(i, j) = r(i, j) + h_der;
      r_minus(i, j) = r(i, j) - h_der;
      f_minus = trialFunction(r_minus, alpha, beta, omega);
      f_plus = trialFunction(r_plus, alpha, beta, omega);
      T_Local = T_Local - (f_minus + f_plus - 2 * f);
      r_plus(i, j) = r(i, j);
      r_minus(i, j) = r(i, j);
    }
  }
  T_Local = 0.5 * h2_der * T_Local / (f);
  return T_Local;
}
