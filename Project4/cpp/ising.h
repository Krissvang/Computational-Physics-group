#ifndef ISING_H
#define	ISING_H

#include <vector>
#include <armadillo>
using namespace std;

void ising(int, int, double, string, double&, double&, double&, double&, double&, vector<int>&, int, int&);
// Function to initialize energy and magnetization
void initialize(int, string, double, arma::Mat<int>&, double&, double&, vector<int>&);
// The metropolis algorithm
void Metropolis(int, arma::Mat<int>&, double&, double&, double *, int&);
#endif /* ISING_H */
