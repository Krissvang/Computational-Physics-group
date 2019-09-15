#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <chrono>
#include <armadillo>

#include "jacobi.h"
// use namespace and armadillo for output and input
using namespace std;
using namespace arma;

int main(int argc, char** argv){
  int n=3;
  double pi=3.14159265358979323846;
  mat A(n,n);
  mat V(n,n);
  double h = 1.0;
  vec r(n);
  int interact = 0;
  double wr = 0;
  vec e_vals(n);
  vec eigenval;


  initialize_beam(n, h, A, r, V, interact, wr);
  eigenval = eig_sym(A);
  jacobi(n,0.001,A,V);
  e_vals=get_eigenvals(A,n);


  cout << e_vals << endl;
  cout << eigenval << endl;
  return 0;
}
