#include <cmath>
#include <iostream>
#include <chrono>
#include <random>

using namespace std;
using namespace std::chrono;

void mc_bruteforce(double (*func) (double *), int n, double R, double& int_mc, double& std_dev, double& time, long t2){
  double x[6], y, fx;
  time_point<high_resolution_clock> start, end;

  start = high_resolution_clock::now();
  int_mc = 0.;  double variance = 0.;
  double sum_f2= 0. ; long idum=t2 ;
  double jacobi_det=pow((2*R),6);
  //Initialize random number generator using mt19937
  random_device rd;
  mt19937_64 gen(rd());
  gen.seed(t2);
  uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
  for ( int i = 1;  i <= n; i++){
//   x[] contains the random numbers for all dimensions
    for (int j = 0; j< 6; j++) {
        x[j]=-R+2*R*RandomNumberGenerator(gen);
    }
    fx=func(x);
    int_mc += fx;
    sum_f2 += fx*fx;
  }
  int_mc = int_mc/((double) n );
  sum_f2 = sum_f2/((double) n );
  variance=sum_f2-int_mc*int_mc;
  end = high_resolution_clock::now();
  duration<double> elapsed = end-start;
  time = elapsed.count();
  int_mc *= jacobi_det;
  std_dev = jacobi_det*sqrt(variance/n);
}

void mc_improved(double (*func) (double *), int n, double& int_mc, double& std_dev, double& time, double& sum_f2, long t2){
  double x[6], y1, y2, r, fx; //x = [r1,r2,theta1,theta2,phi1,phi2]
  time_point<high_resolution_clock> start, end;

  start = high_resolution_clock::now();
  int_mc = 0.;  double variance = 0.;
  sum_f2= 0. ; long idum=t2;
  double jacobi_det = 4*pow(acos(-1.),4.)*1/16;
  //Initialize random number generator using mt19937
  random_device rd;
  mt19937_64 gen(rd());
  gen.seed(t2);
  uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
//   evaluate the integral with importance sampling
  for ( int i = 1;  i <= n; i++){
//   x[] contains the random numbers for all dimensions
    //Generate r1 and r2 according to exponential dist.
    y1 = RandomNumberGenerator(gen);
    y2 = RandomNumberGenerator(gen);
    x[0] = -log(1-y1)/4.;   //-log(1-y1)/(2*alpha)
    x[1] = -log(1-y2)/4.;   //-log(1-y2)/(2*alpha)
    for (int j = 2; j< 4; j++) {
        x[j]=acos(-1.)*RandomNumberGenerator(gen);   //pi*random
    }
    for (int j = 4; j< 6; j++) {
        x[j]=2*acos(-1.)*RandomNumberGenerator(gen); //2*pi*random
    }
    fx=func(x);
    int_mc += fx;
    sum_f2 += fx*fx;
  }
  int_mc = int_mc/((double) n );
  sum_f2 = sum_f2/((double) n );
  variance=sum_f2-int_mc*int_mc;
  end = high_resolution_clock::now();
  duration<double> elapsed = end-start;
  time = elapsed.count();
  int_mc *= jacobi_det;
  std_dev = jacobi_det*sqrt(variance/n);
//   final output
}
