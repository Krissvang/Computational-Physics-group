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
    int n_spins=2; int mcs=1000000;
    double T=1.0;
    double E_avg, heatcap, M_avg, M_abs_avg, susc;
    ising(n_spins,mcs,T,"r",E_avg,heatcap,M_avg,M_abs_avg,susc);
    double Z=12+2*exp(8/T)+2*exp(-8/T);
    REQUIRE(E_avg == Approx((16*exp(-8/T)-16*exp(8/T))/Z).epsilon(0.05));
    REQUIRE(heatcap == Approx(1/(T*T)*(64*(1+3*cosh(8/T)))/(3+cosh(8/T))/(3+cosh(8/T))).epsilon(0.08));
    REQUIRE(susc == Approx(1/T*(32*exp(8/T)+32)/Z).epsilon(0.05));
    REQUIRE(M_abs_avg == Approx((8*exp(8/T)+16)/Z).epsilon(0.05));
    REQUIRE(M_avg == Approx(0.0).epsilon(0.4));
}

TEST_CASE("Testing the metropolis algorithm"){
    
}

