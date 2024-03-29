/*
  This file contains functions to run the variational Monte Carlo algorithm for
  a two-particle system, including functions for local energies and trial
  wavefunctions.
*/

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <armadillo>
#include <iomanip>

using namespace std;
using namespace arma;

//Initialize random number generator using mt19937 algorithm
random_device rd;
mt19937_64 gen(rd());
uniform_real_distribution<double> pos_dist(-1.0, 1.0);
uniform_real_distribution<double> unif_dist(0, 1);
uniform_real_distribution<double> move_dist(-0.5, 0.5);

double Kinetic_E(mat &r, double wfold, double h, double alpha, double beta,
                 double omega, int mcs,
                 double (*trialFunction)(mat &r, double, double, double));

double Potential_E_wo_C(double omega, mat &r);
double Potential_E_w_C(double omega, mat &r);

//Function to initialize positions randomly
mat init_pos()
{
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

//Function to calculate r12=abs(r1-r2)
double r_12(mat &r)
{
  double r_12 = 0;
  for (int i = 0; i < 3; i++)
  {
    r_12 += (r(0, i) - r(1, i)) * (r(0, i) - r(1, i));
  }
  return sqrt(r_12);
}

//Function to execute main part of Monte Carlo simulation
void solver(double &E_avg, double &var_E, double &r12_avg, int &accepted_moves,
            int mcs, double alpha, double beta, double omega,
            double (*TrialWaveFunction)(mat &, double, double, double),
            double (*localE)(mat &, double, double, double), double h,
            double &KE_avg, double &var_KE, double &PE_avg_wo_C, double &var_PE1,
            double &PE_avg_w_C, double &var_PE2)
{
  //Initialize:
  E_avg = r12_avg = KE_avg = var_E = var_KE = PE_avg_wo_C = var_PE1 = PE_avg_w_C = var_PE2 = 0;
  double E2_avg, KE2_avg, PE_w_C_2, PE_wo_C_2 = 0;
  accepted_moves = 0;
  mat r = init_pos();
  double wavefunc = TrialWaveFunction(r, alpha, beta, omega);
  double new_wavefunc;
  double w;
  mat local_r(2, 3);
  double local_energy, local_KE, local_PE_w_C, local_PE_wo_C;
  double local_r12;
  for (int n = 0; n < mcs; n++)
  {
    //Find random new position:
    for (int i = 0; i < 2; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        local_r(i, j) = r(i, j) + h * move_dist(gen);
      }
    }
    //Metropolis algorithm; determine whether move is accepted
    new_wavefunc = TrialWaveFunction(local_r, alpha, beta, omega);
    w = (new_wavefunc * new_wavefunc) / (wavefunc * wavefunc);
    if (unif_dist(gen) <= w)
    {
      r = local_r;
      wavefunc = new_wavefunc;
      accepted_moves += 1;
    }
    //Update energies, distance etc.
    local_KE = Kinetic_E(r, wavefunc, h, alpha, beta, omega, mcs, TrialWaveFunction);
    local_PE_wo_C = Potential_E_wo_C(omega, r);
    local_PE_w_C = Potential_E_w_C(omega, r);
    local_energy = localE(r, alpha, beta, omega);
    local_r12 = r_12(r);

    E_avg += local_energy;
    KE_avg += local_KE;
    PE_avg_w_C += local_PE_w_C;
    PE_avg_wo_C += local_PE_wo_C;

    E2_avg += local_energy * local_energy;
    KE2_avg += local_KE * local_KE;
    PE_w_C_2 += local_PE_w_C * local_PE_w_C;
    PE_wo_C_2 += local_PE_wo_C * local_PE_wo_C;

    r12_avg += local_r12;
  }
  //Calculate averages and variances
  E_avg /= mcs;
  E2_avg /= mcs;
  KE_avg /= mcs;
  KE2_avg /= mcs;
  PE_avg_w_C /= mcs;
  PE_avg_wo_C /= mcs;
  PE_w_C_2 /= mcs;
  PE_wo_C_2 /= mcs;
  r12_avg /= mcs;

  var_E = E2_avg - E_avg * E_avg;
  var_KE = KE2_avg - KE_avg * KE_avg;
  var_PE1 = PE_wo_C_2 - PE_avg_wo_C * PE_avg_wo_C;
  var_PE2 = PE_w_C_2 - PE_avg_w_C * PE_avg_w_C;
}

//Function to find optimal step length h
double FindOptimal_h(double alpha, double beta, double omega, int mcs,
                     double (*TrialWaveFunction)(mat &, double, double,
                                                 double),
                     double (*localE)(mat &, double, double, double))
{
  double h = 3 + .1 / omega;
  double dh = 0.01;
  double error = 0.01;    //1 percent
  int accepted_moves = 0; //To store number of accepted transitions
  double acceptance_ratio = accepted_moves / ((double)mcs);
  bool condition = (0.5 - error) < acceptance_ratio &&
                   acceptance_ratio < (0.5 + error);
  double E_avg, var_E, r12_avg, KE_avg, var_KE, PE1, PE2, Var_PE1, Var_PE2;
  double step = 1;
  while (condition != 1)
  {
    if (h < 0)
    {
      //Reset h and reduce increment
      h = 3 + .1 / omega + step;
      dh *= 0.1;
      step += 1;
    }
    h -= dh;
    solver(E_avg, var_E, r12_avg, accepted_moves, mcs, alpha, beta, omega,
           TrialWaveFunction, localE, h, KE_avg,
           var_KE, PE1, Var_PE1, PE2, Var_PE2);
    acceptance_ratio = accepted_moves / ((double)mcs);
    condition = acceptance_ratio > 0.5 - error &&
                acceptance_ratio < 0.5 + error;
  }
  return h;
}

//Main function, to be called from other programs
void var_mc(double &E_avg, double &var_E, double &r12_avg, int &accepted_moves,
            int mcs, double alpha, double beta, double omega,
            double (*TrialWaveFunction)(mat &, double, double, double),
            double (*localE)(mat &, double, double, double),
            double &KE_avg, double &var_KE,
            double &PE_avg_wo_C, double &var_PE1,
            double &PE_avg_w_C, double &var_PE2)
{
  double h = FindOptimal_h(alpha, beta, omega, 1000, TrialWaveFunction,
                           localE);
  solver(E_avg, var_E, r12_avg, accepted_moves, mcs, alpha, beta, omega,
         TrialWaveFunction, localE, h, KE_avg, var_KE, PE_avg_wo_C, var_PE1,
         PE_avg_w_C, var_PE2);
}

//Function to calculate ri^2 (i= 1 or 2)
double radi(int i, mat &r)
{
  double ri = 0;
  for (int j = 0; j < 3; j++)
  {
    ri += r(i, j) * r(i, j);
  }
  return sqrt(ri);
}

//Function to calculate r1^2+r2^2
double r_squared(mat &r)
{
  double r2 = 0;
  for (int i = 0; i < 3; i++)
  {
    r2 += r(0, i) * r(0, i) + r(1, i) * r(1, i);
  }
  return r2;
}

//Simple trial wavefunction
double TrialWaveFunction1(mat &r, double alpha, double beta, double omega)
{
  double wavefunction = exp(-0.5 * alpha * omega * r_squared(r));
  return wavefunction;
}

//Improved trial wavefunction
double TrialWaveFunction2(mat &r, double alpha, double beta, double omega)
{
  double arg = 1.0 + beta * r_12(r);
  double part2 = exp(r_12(r) / (2 * arg));
  double part1 = TrialWaveFunction1(r, alpha, beta, omega);
  return part1 * part2;
}

//Local energy with "simple" trial wavefunction, no interaction
double E1(mat &r, double alpha, double beta, double omega)
{
  double E = 0.5 * omega * omega * (r_squared(r)) *
                 (1 - alpha * alpha) +
             3 * alpha * omega;
  return E;
}

//Local energy with "simple" trial wavefunction, with interaction
double E_repuls(mat &r, double alpha, double beta, double omega)
{
  double E_repuls = 1 / (r_12(r));
  return E1(r, alpha, beta, omega) + E_repuls;
}

//Local energy for improved trial wavefunction
double E2(mat &r, double alpha, double beta, double omega)
{
  double E = E_repuls(r, alpha, beta, omega) +
             1 / (2 * pow(1 + beta * r_12(r), 2)) *
                 (alpha * omega * r_12(r) - 1 / (2 * pow(1 + beta * r_12(r), 2)) -
                  2 / r_12(r) + 2 * beta / (1 + beta * r_12(r)));
  return E;
}

//Function to calculate kinetic energy
double Kinetic_E(mat &r, double wfold, double h, double alpha, double beta,
                 double omega, int mcs,
                 double (*trialFunction)(mat &r, double, double, double))
{
  double f_minus, f_plus, T_Local;
  mat r_plus(2, 3), r_minus(2, 3);
  r_plus = r_minus = r;
  double h_der = 1. / (double)mcs;   //Step-size for the differentiation
  double h2_der = (double)mcs * mcs; //Inverse squared Step-size for the differentiation
  T_Local = 0;
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      r_plus(i, j) = r(i, j) + h_der;
      r_minus(i, j) = r(i, j) - h_der;
      f_minus = trialFunction(r_minus, alpha, beta, omega);
      f_plus = trialFunction(r_plus, alpha, beta, omega);
      T_Local -= (f_minus + f_plus - 2 * wfold);
      r_plus(i, j) = r(i, j);
      r_minus(i, j) = r(i, j);
    }
  }
  T_Local = 0.5 * h2_der * T_Local / (wfold);

  return T_Local;
}
//Potential energy with out colomb repulsion
double Potential_E_wo_C(double omega, mat &r)
{
  return 0.5 * omega * omega * r_squared(r);

} //Potential energy with colomb repulsion
double Potential_E_w_C(double omega, mat &r)
{
  return 0.5 * omega * omega * r_squared(r) + 1 / r_12(r);
}
