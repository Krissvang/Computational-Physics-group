#include <cmath>
#include <iostream>
#include <chrono>
#include "lib.h"

using namespace std;
using namespace std::chrono;

void mc_bruteforce(double (*func) (double *), int n, double R, double& int_mc, double& std_dev, double& time){
  double x[6], y, fx;
  time_point<high_resolution_clock> start, end;
  time_point<system_clock> time2;
  time2 = system_clock::now();
  start = high_resolution_clock::now();
  duration<double> duration_in_seconds =duration<double>(time2.time_since_epoch());
  long t2= duration_in_seconds.count();
  int_mc = 0.;  double variance = 0.;
  double sum_sigma= 0. ; long idum=t2 ;
  double jacobi_det=pow((2*R),6);

  for ( int i = 1;  i <= n; i++){
//   x[] contains the random numbers for all dimensions
    for (int j = 0; j< 6; j++) {
        x[j]=-R+2*R*ran0(&idum);
    }
    fx=func(x);
    int_mc += fx;
    sum_sigma += fx*fx;
  }
  int_mc = int_mc/((double) n );
  sum_sigma = sum_sigma/((double) n );
  variance=sum_sigma-int_mc*int_mc;
  end = high_resolution_clock::now();
  duration<double> elapsed = end-start;
  time = elapsed.count();
  int_mc *= jacobi_det;
  std_dev = jacobi_det*sqrt(variance/n);
}

void mc_improved(double (*func) (double *), int n, double& int_mc, double& std_dev, double& time, double& sum_sigma, long t2){
  double x[6], y1, y2, r, fx; //x = [r1,r2,theta1,theta2,phi1,phi2]
  time_point<high_resolution_clock> start, end;

  start = high_resolution_clock::now();
  int_mc = 0.;  double variance = 0.;
  sum_sigma= 0. ; long idum=t2;
  double jacobi_det = 4*pow(acos(-1.),4.)*1/16;

//   evaluate the integral with importance sampling
  for ( int i = 1;  i <= n; i++){
//   x[] contains the random numbers for all dimensions
    //Generate r1 and r2 according to exponential dist.
    y1 = ran0(&idum);
    y2 = ran0(&idum);
    x[0] = -log(1-y1)/4.;   //-log(1-y1)/(2*alpha)
    x[1] = -log(1-y2)/4.;   //-log(1-y2)/(2*alpha)
    for (int j = 2; j< 4; j++) {
        x[j]=acos(-1.)*ran0(&idum);   //pi*random
    }
    for (int j = 4; j< 6; j++) {
        x[j]=2*acos(-1.)*ran0(&idum); //2*pi*random
    }
    fx=func(x);
    int_mc += fx;
    sum_sigma += fx*fx;
  }
  int_mc = int_mc/((double) n );
  sum_sigma = sum_sigma/((double) n );
  variance=sum_sigma-int_mc*int_mc;
  end = high_resolution_clock::now();
  duration<double> elapsed = end-start;
  time = elapsed.count();
  int_mc *= jacobi_det;
  std_dev = jacobi_det*sqrt(variance/n);
//   final output
  
}
