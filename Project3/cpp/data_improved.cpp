#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "lib.h"
#include "montecarlo.h"
#include <string>
#include <chrono>

using namespace std::chrono;
using namespace std;

double  improved_MC(double *x);

int main(int argc, char const *argv[]) {
  ofstream ofile;
  int n;
  string filename = "improved_data_";
  double R;
  cout << "Read in the number of Monte-Carlo samples" << endl;
  cin >> n;
  filename.append(to_string(n)+".txt");
  ofile.open(filename);

  n=pow(10,n);
  double int_mc, std_dev, time, sum_sigma;

  time_point<system_clock> time2;
  time2 = system_clock::now();
  duration<double> duration_in_seconds =duration<double>(time2.time_since_epoch());
  long t2= duration_in_seconds.count();

  for (int i = 0; i < atoi(argv[1]); i++)
  {
    mc_improved(&improved_MC,n,int_mc,std_dev,time,sum_sigma,t2+i);
    // final output
    ofile << setprecision(10)<< std_dev <<  setw(30) << setprecision(20) << int_mc << setprecision(10) << setw(20)<< time << endl;
  }
  
  
  ofile.close();
  return 0;
}

double  improved_MC(double *x)
{
   double cosb = cos(x[2])*cos(x[3])+sin(x[2])*sin(x[3])*cos(x[4]-x[5]);
   double deno = sqrt(x[0]*x[0]+x[1]*x[1]-2*x[0]*x[1]*cosb);
   double value = x[0]*x[0]*x[1]*x[1]*sin(x[2])*sin(x[3])/deno;
   return value;

} // end function for the integrand
