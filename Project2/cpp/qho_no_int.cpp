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
using namespace std::chrono;

int main(int argc, char** argv){
  //Checks if enough inputs
    if(argc<=2){
    cout << "Please choose number of mesh points and a max radius."<<endl;
    return 1;
  }
  else{
    //Constants
    int n=atoi(argv[1]);
    double pi=3.14159265358979323846;
    double tol=1e-10;

    //Matrices and vectors needed
    mat A(n,n); mat A2(n,n); mat V(n,n); vec r(n); 
    vec jacobi_e_vals(n); mat e_vecs(3,n);
    double max = atof(argv[2]);
    double h = max/(n+1); double wr = 0;
    int interact = 0;
    vec eigenval;
    mat eigenvec;
    double t_jacobi=0;
    double *t_arma = new double(0);
    int count = 0;
    //Filename
    string filename = "qho_no_int_res_";
    filename.append(to_string(n)+".txt");


    //Initializes the matrices
    initialize_schrodinger(n, h, A, r, V, interact, wr);
    initialize_schrodinger(n, h, A2, r, V, interact, wr);
    
    time_point<high_resolution_clock> start, end;
    start = high_resolution_clock::now();
    //Runs armadillo
    eig_sym(eigenval, eigenvec, A2);
    end = high_resolution_clock::now();
    duration<double> elapsed = end-start;
    *t_arma = elapsed.count();
    //Runs Jacobi algorithm
    jacobi(n,1000,tol,A,V,t_jacobi, count);  
    //Gets the eigenpairs
    jacobi_e_vals = get_eigenvals(A,n);
    e_vecs = get_eigenvecs(A,V,n);
    //Writes to file
    write_eigenpairs(jacobi_e_vals, e_vecs,   eigenval, t_jacobi, *t_arma, filename, count, n);
    return 0;
  }
}