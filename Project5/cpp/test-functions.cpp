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
  double energy, variance, r12;
  var_mc(energy, variance, r12,
         accepted_moves, mcs, alpha, beta, omega,
         TrialWaveFunction1, E1);
  REQUIRE(energy == Approx(3));
  REQUIRE(variance == Approx(0));
  REQUIRE(accepted_moves/((double) mcs) == Approx(0.5).epsilon(0.05));
}

TEST_CASE("Testing functions for r12 and r1^2+r2^2"){
  mat r1 = zeros<mat>(2,3);
  r1(0,0) = 3;
  r1(1,0) = 4;
  mat r2 = zeros<mat>(2,3);
  r2(0,0) = 3;
  r2(1,1) = 4;

  double rsq = r_squared(r1);
  double r12 = r_12(r2);
  REQUIRE(rsq == Approx(25));
  REQUIRE(r12 == Approx(5));
}
