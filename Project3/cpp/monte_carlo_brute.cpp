#include <iostream>
#include <cmath>
#include <chrono>

using namespace std;
using namespace std::chrono;

//     Here we define various functions called by the main program  
//     this function defines the function to integrate  

double func(double x);

//     Main function begins here     
int main()
{
     int n;
     double MCint, MCintsqr2, fx, Variance; 
     cout << "Read in the number of Monte-Carlo samples" << endl;
     cin >> n;
     time_point<high_resolution_clock> start, end;
     start = high_resolution_clock::now();
     MCint = MCintsqr2=0.;
     double invers_period = 1./RAND_MAX; // initialise the random number generator
     srand(time(NULL));  // This produces the so-called seed in MC jargon
//   evaluate the integral with the a crude Monte-Carlo method    
     for ( int i = 1;  i <= n; i++){
  // obtain a floating number x in [0,1]
           double x = double(rand())*invers_period; 
           fx = func(x);
           MCint += fx;
           MCintsqr2 += fx*fx;
     }
     MCint = MCint/((double) n );
     MCintsqr2 = MCintsqr2/((double) n );
     double variance=MCintsqr2-MCint*MCint;
     end = high_resolution_clock::now();
     duration<double> elapsed = end-start;
     double time = elapsed.count();
//   final output 
     cout << "Variance= " << variance << ", Integral = " << MCint << " Exact= " << M_PI  << "Time =" << time << endl;
}  // end of main program 
// this function defines the function to integrate 
double func(double x)
{
  double value;
  value = 4/(1.+x*x);
  return value;
} // end of function to evaluate 