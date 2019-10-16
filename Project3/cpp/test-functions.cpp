#include <cmath>
#include <iostream>
#include <iomanip>
#include "catch.hpp"
#include "lib.h"
#include "montecarlo.h"
#include <chrono>
#define EPS 3.0e-14
#define MAXIT 10

using namespace std::chrono;
using namespace std;

double testfunc_MC(double*);
double testfunc_Gauss(double,double,double,double,double,double);
void gauss_laguerre(double*,double*,int,double);
double gammln(double);


TEST_CASE("Testing improved MC"){
    int n=pow(10,8);
    double int_mc, std_dev, time, sum_sigma;
    time_point<system_clock> time2;
    time2 = system_clock::now();
    duration<double> duration_in_seconds =duration<double>(time2.time_since_epoch());
    long t2= duration_in_seconds.count();
    mc_improved(&testfunc_MC,n,int_mc,std_dev,time,sum_sigma,t2);
    REQUIRE(int_mc==Approx(acos(-1)*acos(-1)/pow(2,6)).epsilon(5*std_dev));
}

TEST_CASE("Testing Gauss-Laguerre"){
  int N = 20;
  double *x_theta = new double [N];
  double *x_phi = new double [N];
  double *x_r = new double [N+1];
  double *w_theta = new double [N];
  double *w_phi = new double [N];
  double *w_r = new double [N+1];

  //Definition of the integration limits for the angular part
  double a_theta = 0;
  double b_theta = 3.14159265359;
  double a_phi = 0;
  double b_phi = 2*3.14159265359;

  //set up the mesh points and weights for the angular part
  gauleg(a_theta,b_theta,x_theta,w_theta, N);
  gauleg(a_phi,b_phi,x_phi,w_phi, N);

  //   set up the mesh points and weights and the power of x^alf
  double alf = 2.0;
  gauss_laguerre(x_r,w_r, N, alf);

  //   evaluate the integral with the Gauss-Laguerre method
  //   Note that we initialize the sum
  double int_gausslag = 0.;
  //six-double loops
  for (int i=1;i<N+1;i++){
      for (int j = 0;j<N;j++){
          for (int k = 0;k<N;k++){
              for (int l = 1;l<N+1;l++){
                  for (int m = 0;m<N;m++){
                      for (int n = 0;n<N;n++){
                          int_gausslag+=w_r[i]*w_theta[j]*w_phi[k]*w_r[l]*w_theta[m]*w_phi[n]*testfunc_Gauss(x_r[i],x_theta[j],x_phi[k],x_r[l],x_theta[m],x_phi[n]);
              }}}}}
  }
  REQUIRE(int_gausslag == Approx(acos(-1)*acos(-1)/pow(2,6)));

  delete [] x_theta;
  delete [] x_phi;
  delete [] x_r;
  delete [] w_theta;
  delete [] w_phi;
  delete [] w_r;
}

double testfunc_MC(double *x){
  return x[0]*x[0]*x[1]*x[1]*sin(x[2])*sin(x[3]);
}

double testfunc_Gauss(double r1, double theta1, double phi1, double r2, double theta2, double phi2)
{
  return sin(theta1)*sin(theta2)/4096.0;
}

void gauss_laguerre(double *x, double *w, int n, double alf)
{
    int i,its,j;
    double ai;
    double p1,p2,p3,pp,z,z1;

    for (i=1;i<=n;i++) {
        if (i == 1) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 2) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-2;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=1;its<=MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    }
}
// end function gaulag

double gammln( double xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

// end function gammln
