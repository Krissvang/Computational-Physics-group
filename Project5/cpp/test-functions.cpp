#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "catch.hpp"
#include "varmontecarlo.h"
#include <armadillo>
using namespace std;

inline int periodic(int i, int limit, int add)
{
  return (i + limit + add) % (limit);
}

TEST_CASE("Testing unperturbed system (energy plus rate of accepted moves)"){
  int mcs = 100000;
  double alpha = 1;
  double omega = 1;
  double beta = 0;
  int accepted_moves;
  double energy, variance, r12, KE, var_KE, PE_wo_C, var_PE_wo, PE_w_C, var_PE_w;
  var_mc(energy, variance, r12,accepted_moves, mcs, alpha, beta, omega,
         TrialWaveFunction1, E1, KE, var_KE, PE_wo_C, var_PE_wo, PE_w_C,
         var_PE_w);
  REQUIRE(energy == Approx(3));
  REQUIRE(variance == Approx(0));
  REQUIRE(accepted_moves/((double) mcs) == Approx(0.5).margin(0.05));
}

TEST_CASE("Testing functions for r12 and r1^2+r2^2"){
  mat r1 = zeros<mat>(2,3);
  mat r2 = zeros<mat>(2,3);
  r1(0,0) = 3;
  r1(1,2) = 4;
  r2(0,0) = 0.1; r2(0,1) = -0.2; r2(0,2) = 0.3;
  r2(1,0) = -0.4; r2(1,1) = -0.5; r2(1,2) = 0.6;

  double rsq1 = r_squared(r1);
  double r121 = r_12(r1);
  double rsq2 = r_squared(r2);
  double r122 = r_12(r2);

  REQUIRE(rsq1 == Approx(25));
  REQUIRE(r121 == Approx(5));
  REQUIRE(rsq2 == Approx(0.91));
  REQUIRE(r122 == Approx(sqrt(0.43)));
}
