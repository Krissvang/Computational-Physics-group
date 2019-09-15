#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <chrono>
#include <armadillo>
#include "jacobi.h"
using namespace std;
using namespace std::chrono;
using namespace arma;

ofstream ofile;

int main(int argc, char* argv[]){
  string fname;         //beginning of filename; filenames will be "fname_m"
  double tol;           //max. allowed off-diagonal value in Jabobi's method
  string tolstr;
  if(argc <= 2){
    cout << "No filename and/or tolerance; read filename and tolerance for "
    "Jacobi method on the same line." << endl;
    exit(1);
  }
  else{
    fname = argv[1];
    tol = atof(argv[2]);
    tolstr = argv[2];
    tolstr.erase(0,2);  //Removes "0." from tol, to avoid period in filename
  }
  int n;
  cout << "Give next value of n (0 to end):" << endl;
  cin >> n;

  while(n != 0){
    string outfile = fname;
    outfile.append("_"+to_string(n)+"_"+tolstr);

    mat A(n,n);
    mat V(n,n);
    double h = 1./n;
    vec r(n);
    int interact = 0;
    double wr = 0;
    vec e_vals(n);
    vec eigenval;
    double t_Jacobi = 0;
    double t_arma = 0;
    time_point<high_resolution_clock> start_arma, finish_arma;

    initialize_beam(n, h, A, r, V, interact, wr);
    start_arma = high_resolution_clock::now();
    eigenval = eig_sym(A);
    finish_arma = high_resolution_clock::now();
    duration<double> elapsed = finish_arma-start_arma;
    t_arma = elapsed.count();
    jacobi(n,tol,A,V,t_Jacobi);
    e_vals=get_eigenvals(A,n);

    ofile.open(outfile);
    ofile << setiosflags(ios::showpoint);
    ofile << "t_Jacobi: " << scientific << t_Jacobi << " t_armadillo: "
    << scientific << t_arma << endl;
    ofile << "Jacobi eigenvalues:    Armadillo eigenvalues:" << endl;
    for(int i = 0; i < n; i++){
      ofile << setw(19) << setprecision(8) << e_vals[i];
      ofile << setw(26) << setprecision(8) << eigenval[i] << endl;
    }
    ofile.close();

    cout << "Give next value of n (0 to end): " << endl;
    cin >> n;
  }
  return 0;
}
