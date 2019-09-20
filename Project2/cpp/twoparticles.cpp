#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <chrono>
#include <armadillo>
#include "jacobi.h"
using namespace std;
using namespace arma;

ofstream ofile;

int main(int argc, char const *argv[]){
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
    fname = "results/"+fname;
    tol = atof(argv[2]);
    tolstr = argv[2];
  }
  int n;
  cout << "Give next value of n (0 to end):" << endl;
  cin >> n;
  while(n != 0){
    string outfile = fname;
    outfile.append("_"+to_string(n)+"_"+tolstr);

    mat A(n,n);
    mat V(n,n);
    vec r(n);
    int interact = 1;
    double wr = 0;
    vec e_vals(n);
    mat e_vecs(3,n);
    double t = 0;
    int count;
    vec wrlist = vec("0.01 0.5 1.0 5.0");
    vec max_r_list(4);
    for(int i=0; i < 4; i++){
      cout << "Give max. r/alpha for wr = " << wrlist[i] << ": " << endl;
      cin >> max_r_list(i);
    }
    for(int i=0; i < 4; i++){
      double max_r = max_r_list[i];
      double h = max_r/(n+1);
      initialize_schrodinger(n, h, A, r, V, interact, wrlist[i]);
      jacobi(n,tol,A,V,t,count);
      e_vals=get_eigenvals(A,n);
      e_vecs=get_eigenvecs(A,V,n);
      ofile.open(outfile+"_"+to_string(int(max_r))+"_"+to_string(i));
      ofile << setiosflags(ios::showpoint);
      ofile << "wr: " << wrlist[i] << ". time: " << scientific << t
      << " (" << count << " iterations)" << endl;
      ofile << "Jacobi eigenvalues:       Ground state:" << endl;
      for(int i = 0; i < n; i++){
        ofile << setw(19) << setprecision(8) << e_vals[i];
        ofile << setw(19) << setprecision(8) << e_vecs(0,i) << endl;
    }
    ofile.close();
    }
    cout << "Give next value of n (0 to end): " << endl;
    cin >> n;
  }
  return 0;
}
