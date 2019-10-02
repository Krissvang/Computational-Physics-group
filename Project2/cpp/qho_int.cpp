#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <chrono>
#include <armadillo>
#include <sstream>

#include "jacobi.h"
// use namespace and armadillo for output and input
using namespace std;
using namespace arma;
using namespace std::chrono;

int main(int argc, char** argv){
  //Checks if enough inputs
    if(argc<=3){
    cout << "Please choose number of mesh points, a max radius and a frequency"<<endl;
    return 1;
  }
  else{

    //Constants
    int n=atoi(argv[1]);
    double pi=3.14159265358979323846;
    double tol=1e-16;
    

    //Matrices, vectors and variables needed
    mat A(n,n); mat A2(n,n); mat V(n,n); vec r(n);
    mat e_vecs(3,n); vec jacobi_e_vals(n);
    double max = atof(argv[2]);
    double h = max/(n+1); double wr = atof(argv[3]);
    int interact = 1;
    vec eigenval;
    mat eigenvec;
    double t_jacobi=0;
    double *t_arma = new double(0);
    int count = 0;
    stringstream ss;
    ss << wr;
    string filename = "qho_int_res_";
    filename.append(to_string(n)+"_wr="+ss.str()+".txt");
    


    //Initializing matrices
    initialize_schrodinger(n, h, A, r, V, interact, wr);
    initialize_schrodinger(n, h, A2, r, V, interact, wr);
    //Takes the time
    time_point<high_resolution_clock> start, end;
    start = high_resolution_clock::now();
    //run armadillo
    eig_sym(eigenval, eigenvec, A2);
    end = high_resolution_clock::now();
    duration<double> elapsed = end-start;
    *t_arma = elapsed.count();
    //Runs the Jacobi algorithm
    jacobi(n,1000,tol,A,V,t_jacobi, count);  
    //Gets the eigenpairs
    jacobi_e_vals = get_eigenvals(A,n);
    e_vecs = get_eigenvecs(A,V,n);
    //Writes to file
    write_eigenpairs(jacobi_e_vals, e_vecs,   eigenval, t_jacobi, *t_arma, filename , count, n);
    return 0;
  }
}