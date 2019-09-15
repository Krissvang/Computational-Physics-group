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
  //Constants
  int n=3;
  double pi=3.14159265358979323846;

  //Matrices and vectors needed
  mat A(n,n); mat A2(n,n); mat V(n,n); vec r(n);
  double h = 1.0; double wr = 0;
  int interact = 0;
  vec eigenval;
  mat eigenvec;

  //Just testing the algorithms
  initialize_classic(n, h, A, r, V, interact, wr);
  initialize_classic(n, h, A2, r, V, interact, wr);
  cout << V << endl;
  eig_sym(eigenval, eigenvec, A);
  jacobi(n,20,0.0000000000000001,A,V);
  cout << eigenvec << endl;
  cout << V << endl;

  return 0;
}