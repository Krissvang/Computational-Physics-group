/*
   Program to solve the two-dimensional Ising model
   with zero external field.
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
*/

#include <cmath>
#include <random>
#include <armadillo>
using namespace  std;
using namespace arma;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) {
  return (i+limit+add) % (limit);
}

// Function to initialize energy and magnetization
void initialize(int, string, double, arma::Mat<int>&, double&, double&, vector<int>&);
// The metropolis algorithm
void Metropolis(int, arma::Mat<int>&, double&, double&, double *);

//Initialize random number generator using mt19937
random_device rd;
mt19937_64 gen(rd());
uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);

/* function for finding the various properties (avg. energy, magnetization etc.)
   for a given temperature, lattice size, and #Monte Carlo cycles
*/
void ising(int n_spins, int mcs, double temperature, string init, double& E_avg,
double& heatcap, double& M_avg, double& M_abs_avg, double& susc,
vector<int>& count, int steady_start)
{
  Mat<int> spin_matrix(n_spins,n_spins);
  double exp_de[17], E, M, E2_avg, M2_avg;
  double Emin = -2*n_spins*n_spins;
  //    initialise energy and magnetization
  E = M = E_avg = E2_avg = M_avg = M2_avg = M_abs_avg = 0.;
  // setup array for possible energy changes
  for( int de =-8; de <= 8; de++) exp_de[de+8] = 0;
  for( int de =-8; de <= 8; de+=4) exp_de[de+8] = exp(-de/temperature);
  // initialise array for expectation values
  initialize(n_spins, init, temperature, spin_matrix, E, M, count);
  // start Monte Carlo computation
  // First reach equilibrium:
  for (int cycles = 1; cycles <= steady_start; cycles++){
    Metropolis(n_spins, spin_matrix, E, M, exp_de);
  }
  // Then collect data
  for (int cycles = steady_start+1; cycles <= mcs; cycles++){
      Metropolis(n_spins, spin_matrix, E, M, exp_de);
      // update expectation values
      E_avg += E; E2_avg += E*E; M_avg += M;
      M2_avg += M*M; M_abs_avg += fabs(M);
      count[(E-Emin)/4] += 1;
  }
  double mcs_inv = 1./((double) (mcs-steady_start));
  E_avg *= mcs_inv; E2_avg *= mcs_inv; M_avg *= mcs_inv;
  M2_avg *= mcs_inv; M_abs_avg *= mcs_inv;
  double E_var = E2_avg-E_avg*E_avg;
  double M_var = M2_avg-M_avg*M_avg;
  heatcap = E_var/(temperature*temperature);
  susc = M_var/temperature;
}


// function to initialize energy, spin matrix and magnetization, and counter
void initialize(int n_spins, string init, double temperature,
  Mat<int>& spin_matrix, double& E, double& M, vector<int>& count)
{
  // setup spin matrix and intial magnetization
  for(int i =0; i < n_spins; i++) {
    for (int j= 0; j < n_spins; j++){
      if(init == "r"){
        double r = -1+2*RandomNumberGenerator(gen);
        spin_matrix(i,j) = (int) (r/fabs(r));// Generate initial state randomly
      }
      else{
        spin_matrix(i,j) = 1;      //Ordered initial state
      }
      M +=  (double) spin_matrix(i,j);
    }
  }
  // setup initial energy
  for(int i =0; i < n_spins; i++) {
    for (int j= 0; j < n_spins; j++){
      E -=  (double) spin_matrix(i,j)*
	(spin_matrix(periodic(i,n_spins,-1),j) +
	 spin_matrix(i,periodic(j,n_spins,-1)));
    }
  }
  // setup counter
  int length = count.size();
  for(int i = 0; i < length; i++){
    count[i] = 0;
  }
}// end function initialize

void Metropolis(int n_spins, Mat<int>& spin_matrix, double& E, double&M, double *exp_de)
{
  // loop over all spins
  for(int i =0; i < n_spins; i++) {
    for (int j= 0; j < n_spins; j++){
      int irand = (int) (RandomNumberGenerator(gen)*(double)n_spins);
      int jrand = (int) (RandomNumberGenerator(gen)*(double)n_spins);
      int deltaE =  2*spin_matrix(irand,jrand)*
	(spin_matrix(irand,periodic(jrand,n_spins,-1))+
	 spin_matrix(periodic(irand,n_spins,-1),jrand) +
	 spin_matrix(irand,periodic(jrand,n_spins,1)) +
	 spin_matrix(periodic(irand,n_spins,1),jrand));
      if (RandomNumberGenerator(gen) <= exp_de[deltaE+8] ) {
	spin_matrix(irand,jrand) *= -1;  // flip one spin and accept new spin config
        M += (double) 2*spin_matrix(irand,jrand);
        E += (double) deltaE;
      }
    }
  }
} // end of Metropolis sampling over spins
