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
    if(argc<=2){
    cout << "Please choose number of mesh points and a max radius."<<endl;
    return 1;
  }
  else{
    //Constants
    int n=atoi(argv[1]);
    double pi=3.14159265358979323846;
    double tol=1e-16;

    //Matrices and vectors needed
    double max = atof(argv[2]);
    double h = max/(n+1); double wr = 0;
    int interact = 0;
    vec eigenval;
    mat eigenvec;
    double t_jacobi=0;
    int count = 0;

    



    int *p_i = new int (200);

    double *pe1= new double, *pe2=new double, *pe3=new double;
 
    do{
      mat A(*p_i,*p_i); mat V(*p_i,*p_i); vec r(*p_i);
      cout << *p_i << endl;
      initialize_schrodinger(*p_i, h, A, r, V, interact, wr);
      jacobi(*p_i,1000,tol,A,V,t_jacobi, count);  
      vec jacobi_e_vals(*p_i);
      jacobi_e_vals = get_eigenvals(A,*p_i);
      *pe1 = jacobi_e_vals(0);
      *pe2 = jacobi_e_vals(1);
      *pe3 = jacobi_e_vals(2);
      *p_i+=50;
      cout << abs(3.0-*pe1) << endl;
    }
    while(abs(3.0-*pe1)>1e-4);
    string filename = "output/qho_leading_digits.txt";
    ofstream ofile;
    ofile.open(filename);
    ofile << setiosflags(ios::showpoint);
    ofile << "Number of iterations" << endl;
    ofile << count << endl;
    ofile << "Eigenvalues"<< endl;
    ofile << *pe1 << endl;
    ofile << *pe2 << endl;
    ofile << *pe3 << endl;
    ofile.close();
    cout << *pe1 << endl;
    delete pe1, pe2, pe3, p_i;
    return 0;
  }
}