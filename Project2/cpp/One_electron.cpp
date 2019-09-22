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
        //fname = "results/"+fname;
        tol = atof(argv[2]);
        tolstr = argv[2];
    }
    int n;
    cout << "Give next value of n (0 to end):" << endl;
    cin >> n;
    
    while(n != 0){
        string outfile = fname;
        //outfile.append("_"+to_string(n)+"_"+tolstr);
        
        mat A(n,n);
        mat V(n,n);
        double max_r;
        double h;
        vec r(n);
        int interact = 0;
        double wr = 0;
        vec e_vals(n);
        int n_eigenval=5;
        vec eigenval(n_eigenval);
        //mat e_vecs(3,n);
        
        if(n<n_eigenval){
            cout<<"n must be bigger than 4"<<endl;
            exit(1);
        }
        
        cout<<"Give max. r/alpha:"<<endl;
        cin>>max_r;
        h=max_r/(n+1);
        
        initialize(n, h, A, r, V, interact, wr);
        //cout<<A<<endl;
        jacobi(n,tol,A,V);
        //cout<<A<<endl;
        e_vals=get_eigenvals(A,n);
        //e_vecs=get_eigenvecs(A,V,n);
        //cout<<e_vals;
        
        for(int i=0; i<n_eigenval; i++){
            eigenval(i)=e_vals(i);
        }
        
        //cout<<eigenval;
        ofile.open(outfile+".txt");
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << "r_max:" << max_r << endl;
        ofile << "Jacobi eigenvalues:"<<endl;
        for(int i=0; i<n_eigenval; i++){
            ofile << setw(19) << setprecision(8) << eigenval[i] <<endl;
        }
        
        ofile.close();
        
        cout << "Give next value of n (0 to end): " << endl;
        cin >> n;
    }
    
    
    return 0;
}
