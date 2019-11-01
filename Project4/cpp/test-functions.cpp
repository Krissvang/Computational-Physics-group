#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include "catch.hpp"
#include "ising.h"
#include <random>
#include <armadillo>
using namespace std;
using namespace arma;

TEST_CASE("Testing the ising model"){
    int n_spins=2; int mcs=1000;
    double T=1.0;
    double E_avg, heatcap, M_avg, M_abs_avg, susc;
    ising(n_spins,mcs,temperature,init,E_avg,heatcap,M_avg,M_abs_avg,susc);
    double Z=12+2*exp(8/T)+2*exp(-8/T);
    REQUIRE(E_avg==Approx(16*exp(8/T)-16*exp(-8/T)).epsilon(0.1));
}

