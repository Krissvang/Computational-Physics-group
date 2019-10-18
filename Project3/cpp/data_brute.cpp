#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <random>
#include "montecarlo.h"
#include <string>
#include <chrono>

using namespace std::chrono;
using namespace std;

double brute_force_MC(double *);

int main(int argc, char const *argv[]) {
  ofstream ofile;
  int n;
  string filename = "brute_data_";
  double R;
  cout << "Read in the number of Monte-Carlo samples" << endl;
  cin >> n;
  filename.append(to_string(n)+".txt");
  ofile.open(filename);

  n=pow(10,n);
  cout << "Read in R (max. absolute value of x_i)" << endl;
  cin >> R;
  double int_mc, std_dev, time;

  time_point<system_clock> time2;
  time2 = system_clock::now();
  duration<double> duration_in_seconds =duration<double>(time2.time_since_epoch());
  long t2= duration_in_seconds.count();

  for (int i = 0; i < atoi(argv[1]); i++)
  {
    mc_bruteforce(&brute_force_MC,n,R,int_mc,std_dev,time,t2+i);
    // final output
    ofile << setw(20)<< std_dev <<  setw(20)<<setprecision(10) << int_mc << setw(20) <<setprecision(10)<< time << endl;
  }


  ofile.close();
  return 0;
}

// this function defines the integrand to integrate
double  brute_force_MC(double *x)
{
   double alpha = 2.;
// evaluate the different terms of the exponential
   double exp1=-2*alpha*sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
   double exp2=-2*alpha*sqrt(x[3]*x[3]+x[4]*x[4]+x[5]*x[5]);;
   double deno = sqrt(pow((x[0]-x[3]),2)+pow((x[1]-x[4]),2)+pow((x[2]-x[5]),2));
   return exp(exp1+exp2)/deno;
} // end function for the integrand