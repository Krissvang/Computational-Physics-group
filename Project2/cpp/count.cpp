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
    if(argc<=1){
    cout << "Please choose number of mesh points"<<endl;
    return 1;
  }
  else{

    //Constants
    int n=atoi(argv[1]);
    double pi=3.14159265358979323846;
    double tol=1e-16;

    //Matrices and vectors needed
    mat A(n,n); mat A2(n,n); mat V(n,n); vec r(n);
    double h = 1.0/(n+1); double wr = 0;
    int interact = 0;
    vec eigenval;
    mat eigenvec;
    double t_jacobi=0;
    int count = 0;
    output_count_beam(atoi(argv[1]), "number_of_transfos.txt",h,wr,interact);  
    return 0;
  }
}