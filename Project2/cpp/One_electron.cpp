//
//  One_electron.cpp
//  
//
// 
//

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
        double h = 1./(n+1);
        vec r(n);
        int interact = 0;
        double wr = 0;
        vec e_vals(n);
        vec eigenval;
        int count;
        
        initialize(n, h, A, r, V, interact, wr);
        //cout<<A<<endl;
        jacobi(n,tol,A,V,t_Jacobi,count);
        e_vals=get_eigenvals(A,n);
        e_vecs=get_eigenvecs(A,V,n);
        
        cout << "Give next value of n (0 to end): " << endl;
        cin >> n;
    }
    
    
    return 0;
}
