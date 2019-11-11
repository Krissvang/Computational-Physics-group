<<<<<<< HEAD
#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "catch.hpp"
#include "ising.h"
#include <armadillo>
using namespace std;

inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}

TEST_CASE("Testing the ising model"){
    int n_spins=2; int mcs=1000000; int steady_start=20;  int count_configs = 0;
    vector<int> count(n_spins * n_spins + 1);
    double T=1.0;
    double E_avg, heatcap, M_avg, M_abs_avg, susc;
    ising(n_spins,mcs,T,"r",E_avg,heatcap,M_avg,M_abs_avg,susc,count,steady_start,count_configs);
    double Z=12+2*exp(8/T)+2*exp(-8/T);
    REQUIRE(E_avg == Approx((16*exp(-8/T)-16*exp(8/T))/Z).epsilon(0.05));
    REQUIRE(heatcap == Approx(1/(T*T)*(64*(1+3*cosh(8/T)))/(3+cosh(8/T))/(3+cosh(8/T))).epsilon(0.1));
    REQUIRE(susc == Approx(1/T*(32*exp(8/T)+32)/Z).epsilon(0.05));
    REQUIRE(M_abs_avg == Approx((8*exp(8/T)+16)/Z).epsilon(0.05));
    REQUIRE(M_avg+1 == Approx(1).epsilon(0.7));
}

TEST_CASE("Testing the initialize algorithm"){
    int n_spins=2; int count_configs = 0;
    double T=1;
    arma::Mat<int> spin_matrix(n_spins,n_spins);
    vector<int> count;
    double exp_de[17], E, M, E_test;
    E = M = 0.;
    for( int de =-8; de <= 8; de++) exp_de[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) exp_de[de+8] = exp(-de/T);
    initialize(n_spins, "a", T, spin_matrix, E, M, count);
    int sum_elements=spin_matrix(0,0)+spin_matrix(1,0)+spin_matrix(0,1)+spin_matrix(1,1);
    for(int i =0; i < n_spins; i++) {
      for (int j= 0; j < n_spins; j++){
        E_test -=  (double) spin_matrix(i,j)*
      (spin_matrix(periodic(i,n_spins,-1),j) +
       spin_matrix(i,periodic(j,n_spins,-1)));
      }}
    REQUIRE(sum_elements<=4);
    REQUIRE(sum_elements>=-4);
    REQUIRE(E_test==E);
    REQUIRE_FALSE(spin_matrix(0,0)==0);
    REQUIRE_FALSE(spin_matrix(1,0)==0);
    REQUIRE_FALSE(spin_matrix(0,1)==0);
    REQUIRE_FALSE(spin_matrix(1,1)==0);
}


TEST_CASE("Testing the metropolis algorithm"){
    int n_spins=2; int count_configs = 0;
    double T=1;
    arma::Mat<int> spin_matrix(n_spins,n_spins);
    vector<int> count;
    double exp_de[17], E, M;
    E = M = 0.;
    for( int de =-8; de <= 8; de++) exp_de[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) exp_de[de+8] = exp(-de/T);
    initialize(n_spins, "r", T, spin_matrix, E, M, count);
    Metropolis(n_spins, spin_matrix, E, M, exp_de, count_configs);
    int sum_elements=spin_matrix(0,0)+spin_matrix(1,0)+spin_matrix(0,1)+spin_matrix(1,1);
    REQUIRE(count_configs>=0);
    REQUIRE(count_configs<=4);
    REQUIRE(sum_elements<=4);
    REQUIRE(sum_elements>=-4);
    REQUIRE_FALSE(spin_matrix(0,0)==0);
    REQUIRE_FALSE(spin_matrix(1,0)==0);
    REQUIRE_FALSE(spin_matrix(0,1)==0);
    REQUIRE_FALSE(spin_matrix(1,1)==0);
    
}

//=======
#include <cmath>
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "catch.hpp"
#include "ising.h"
#include <armadillo>
using namespace std;

inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}

TEST_CASE("Testing the ising model"){
    int n_spins=2; int mcs=1000000; int steady_start=20;  int count_configs = 0;
    vector<int> count(n_spins * n_spins + 1);
    double T=1.0;
    double E_avg, heatcap, M_avg, M_abs_avg, susc;
    ising(n_spins,mcs,T,"r",E_avg,heatcap,M_avg,M_abs_avg,susc,count,steady_start,count_configs);
    double Z=12+2*exp(8/T)+2*exp(-8/T);
    REQUIRE(E_avg == Approx((16*exp(-8/T)-16*exp(8/T))/Z).epsilon(0.05));
    REQUIRE(heatcap == Approx(1/(T*T)*(64*(1+3*cosh(8/T)))/(3+cosh(8/T))/(3+cosh(8/T))).epsilon(0.1));
    REQUIRE(susc == Approx(1/T*(32*exp(8/T)+32)/Z).epsilon(0.05));
    REQUIRE(M_abs_avg == Approx((8*exp(8/T)+16)/Z).epsilon(0.05));
    REQUIRE(M_avg+1 == Approx(1).epsilon(0.7));
}

TEST_CASE("Testing the initialize algorithm"){
    int n_spins=2; int count_configs = 0;
    double T=1;
    arma::Mat<int> spin_matrix(n_spins,n_spins);
    vector<int> count;
    double exp_de[17], E, M, E_test;
    E = M = 0.;
    for( int de =-8; de <= 8; de++) exp_de[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) exp_de[de+8] = exp(-de/T);
    initialize(n_spins, "a", T, spin_matrix, E, M, count);
    int sum_elements=spin_matrix(0,0)+spin_matrix(1,0)+spin_matrix(0,1)+spin_matrix(1,1);
    for(int i =0; i < n_spins; i++) {
      for (int j= 0; j < n_spins; j++){
        E_test -=  (double) spin_matrix(i,j)*
      (spin_matrix(periodic(i,n_spins,-1),j) +
       spin_matrix(i,periodic(j,n_spins,-1)));
      }}
    REQUIRE(sum_elements<=4);
    REQUIRE(sum_elements>=-4);
    REQUIRE(E_test==E);
}


TEST_CASE("Testing the metropolis algorithm"){
    int n_spins=2; int count_configs = 0;
    double T=1;
    arma::Mat<int> spin_matrix(n_spins,n_spins);
    vector<int> count;
    double exp_de[17], E, M;
    E = M = 0.;
    for( int de =-8; de <= 8; de++) exp_de[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) exp_de[de+8] = exp(-de/T);
    initialize(n_spins, "r", T, spin_matrix, E, M, count);
    Metropolis(n_spins, spin_matrix, E, M, exp_de, count_configs);
    int sum_elements=spin_matrix(0,0)+spin_matrix(1,0)+spin_matrix(0,1)+spin_matrix(1,1);
    REQUIRE(count_configs>=0);
    REQUIRE(count_configs<=4);
    REQUIRE(sum_elements<=4);
    REQUIRE(sum_elements>=-4);
    
}

>>>>>>> be60cecf53c15ecf44424cdec0c8010d5c77939a
