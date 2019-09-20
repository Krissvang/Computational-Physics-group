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
  int n=100;
  double pi=3.14159265358979323846;
  double tol=1e-16;

  //Matrices and vectors needed
  mat A(n,n); mat A2(n,n); mat V(n,n); vec r(n);
  double max = 10;
  double h = max/(n+1); double wr = 0;
  int interact = 0;
  vec eigenval;
  mat eigenvec;
  double t_jacobi=0;
  int count = 0;

  


  //Just testing the algorithms
  initialize_schrodinger(n, h, A, r, V, interact, wr);
  initialize_schrodinger(n, h, A2, r, V, interact, wr);

  eig_sym(eigenval, eigenvec, A2);

  jacobi(n,1000,tol,A,V,t_jacobi, count);  
  vec jacobi_e_vals(n);
  jacobi_e_vals = get_eigenvals(A,n);
  mat e_vecs(3,n);
  e_vecs = get_eigenvecs(A,V,n);

  string filename = "qho_no_int_res_";
  filename.append(to_string(n)+".txt");
  
  write_eigenpairs(jacobi_e_vals, e_vecs,   eigenval, t_jacobi, 0, filename, count, n);

  //output_count_beam(200, "number_of_transfos.txt",h,wr,interact);

  
  return 0;
}