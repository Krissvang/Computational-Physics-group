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
    

    string outfile = fname;
    //outfile.append("_"+to_string(n)+"_"+tolstr);
    
    int selection;
    cout<<"Write 1 to keep n constant and to vary r_max. Write 0 for the defaul program"<<endl;
    cin>>selection;
    
    if(selection==0){
            int n;
            cout << "Give the value of n (0 to end):" << endl;
            cin >> n;
    
            while(n != 0){
        
                mat A(n,n);
                mat V(n,n);
                double max_r;
                double h;
                vec r(n);
                int interact = 0;
                double wr = 0;
                vec e_vals(n);
                mat e_vecs(3,n);
        
                cout<<"Give max. r/alpha:"<<endl;
                cin>>max_r;
                h=max_r/(n+1);
        
                initialize(n, h, A, r, V, interact, wr);
                //cout<<A<<endl;
                jacobi(n,tol,A,V);
                //cout<<A<<endl;
                e_vals=get_eigenvals(A,n);
                e_vecs=get_eigenvecs(A,V,n);
        
                //cout<<eigenval;
                ofile.open(outfile);
                ofile << setiosflags(ios::showpoint | ios::uppercase);
                ofile << "r_max:" << max_r << endl;
                ofile << "tolerance:" << tol << endl;
                ofile << "n:" << n << endl;
                ofile << "Jacobi eigenvalues:           v0:               v1:              v2:"<<endl;
                for(int i=0; i<n; i++){
                    ofile << setw(19) << setprecision(8) << e_vals[i];
                    ofile << setw(19) << setprecision(8) << e_vecs(0,i);
                    ofile << setw(19) << setprecision(8) << e_vecs(1,i);
                    ofile << setw(19) << setprecision(8) << e_vecs(2,i) << endl;
                }
        
                ofile.close();
        
                cout << "Give next value of n (0 to end): " << endl;
                cin >> n;
            }
        }
            
            
    if(selection==1){
            
            int n;
            cout << "Give the value of n (0 to end):" << endl;
            cin >> n;
        
            if(n==0) exit(1);
        
            mat A(n,n);
            mat V(n,n);
            double max_r;
            double min_r;
            double h;
            vec r(n);
            int interact = 0;
            int r_interval=0;
            double wr = 0;
            vec r_rep(10);
            vec e_vals(n);
            mat e_vecs(3,n);
            
            cout<<"Give max. r/alpha:"<<endl;
            cin>>max_r;
            
            cout<<"Give min. r/alpha:"<<endl;
            cin>>min_r;
            
            r_interval=(max_r-min_r)/10;
            r_rep(0)=min_r;
            for(int i=1; i<10; i++){
                r_rep(i)=r_rep(i-1)+r_interval;
            }
            
            for (int j=0; j<10; j++){
                h=r_rep(j)/(n+1);
        
                initialize(n, h, A, r, V, interact, wr);
                jacobi(n,tol,A,V);
                e_vals=get_eigenvals(A,n);
                e_vecs=get_eigenvecs(A,V,n);
                
                string fileout = fname;
                string argument = to_string(j);
                fileout.append(argument);
                
                ofile.open(fileout);
                ofile << setiosflags(ios::showpoint | ios::uppercase);
                ofile << "r:" << r_rep[j] << endl;
                ofile << "tolerance:" << tol << endl;
                ofile << "n:" << n << endl;
                ofile << "Jacobi eigenvalues:           v0:               v1:             v2:"<<endl;
                for(int i=0; i<n; i++){
                    ofile << setw(19) << setprecision(8) << e_vals[i];
                    ofile << setw(19) << setprecision(8) << e_vecs(0,i);
                    ofile << setw(19) << setprecision(8) << e_vecs(1,i);
                    ofile << setw(19) << setprecision(8) << e_vecs(2,i) << endl;
                }
        
            ofile.close();
            }
        }
    
    return 0;
}
