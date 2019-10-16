#include <cmath>
#include <iostream>
#include <iomanip>
#include "catch.hpp"
#include "lib.h"
#include "montecarlo.h"
#include <chrono>

using namespace std::chrono;
using namespace std;

double testfunc_MC(double *x){
  return x[0]*x[0]*x[1]*x[1]*sin(x[2])*sin(x[3]);
}


TEST_CASE("Testing improved MC"){
    int n=pow(10,8);
    double int_mc, std_dev, time, sum_sigma;
    time_point<system_clock> time2;
    time2 = system_clock::now();
    duration<double> duration_in_seconds =duration<double>(time2.time_since_epoch());
    long t2= duration_in_seconds.count();
    mc_improved(&testfunc_MC,n,int_mc,std_dev,time,sum_sigma,t2);
    int R = 4;
    REQUIRE(int_mc==Approx(acos(-1)*acos(-1)/pow(2,6)).epsilon(5*std_dev));
}
