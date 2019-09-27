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

//Begin main program
int main(int argc, char* argv[]){
    string fname;         //beginning of filename; filenames will be "fname_m"
    double tol;           //max. allowed off-diagonal value in Jabobi's method
    string tolstr;
    double time=0;
    int count=0;
    
    //We get the filename and the tolerance from the user
    if(argc <= 2){
        cout << "No filename and/or tolerance; read filename and tolerance for "
        "Jacobi method on the same line." << endl;
        exit(1);
    }
    else{
        fname= "output/";
        fname += argv[1];
        //fname = "results/"+fname;
        tol = atof(argv[2]);
        tolstr = argv[2];
    }
    

    //Definition of the output file
    string outfile = fname;
    
    //We let the user choose if he wants to vary r_max or not
    int selection;
    cout<<"Write 1 to keep n constant and to vary r_max. Write 0 for the defaul program"<<endl;
    cin>>selection;
    
    //Case in which we do not vary n
    if(selection==0){
            int n;
            //Here the user chooses the number of points of the discretize solution
            cout << "Give the value of n (0 to end):" << endl;
            cin >> n;
        
            //While cicle to make the user choose different values of n after the first iteration
            while(n != 0){
            
                //declaration of varibles
                mat A(n,n); mat V(n,n); mat e_vecs(3,n);
                double max_r; double h; double wr = 0;
                vec r(n); vec e_vals(n);
                int interact = 0;
                
                //Here the user chooses the value of r_max
                cout<<"Give max. r/alpha:"<<endl;
                cin>>max_r;
                //Definition of the step lenght
                h=max_r/(n+1);
        
                //Initialization of the matrix (this function are written in the jacobi.cpp file)
                initialize_schrodinger(n, h, A, r, V, interact, wr);
                
                //Diagonalization of the matrix (this function are written in the jacobi.cpp file)
                jacobi(n,1000,tol,A,V,time,count);
                //cout<<A<<endl;
                
                //Here we get the eigenvalues and eigenvectors (this function are written in the jacobi.cpp file)
                e_vals=get_eigenvals(A,n);
                e_vecs=get_eigenvecs(A,V,n);
        
                //Here we print the results on a file
                //We open the file
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
        
                //We close the file
                ofile.close();
        
                //Here we let the user choose if he wants to run the program with a different value for n
                cout << "Give next value of n (0 to end): " << endl;
                cin >> n;
            }
        }
           
    //Case in which we vary n
    if(selection==1){
            
            int n;
            //Here the user chooses the number of points of the discretize solution
            cout << "Give the value of n (0 to end):" << endl;
            cin >> n;
        
            if(n==0) exit(1);
        
            //declaration of varibles
            mat A(n,n); mat V(n,n); mat e_vecs(3,n);
            double max_r; double min_r; double h; double r_interval=0; double wr = 0;
            vec r(n); vec r_rep(10); vec e_vals(n);
            int interact = 0;
            
            //Here the user chooses the interval in which to vary r_max
            cout<<"Give max. r/alpha:"<<endl;
            cin>>max_r;
            
            cout<<"Give min. r/alpha:"<<endl;
            cin>>min_r;
            
            //Definition of the step between two different values of r_max
            r_interval=(max_r-min_r)/double(10.0);
            r_rep(0)=min_r;
            //Definition of the vector for the r_max values. For each therm of the vector we will do an iteration
            for(int i=1; i<10; i++){
                r_rep(i)=r_rep(i-1)+r_interval;
            }
            
            //Cicle on the the different values of r_max (we have defined 10 possible values fot this variable)
            for (int j=0; j<10; j++){
                //Definition of the step lenght
                h=r_rep(j)/(n+1);
        
                initialize_schrodinger(n, h, A, r, V, interact, wr);
                jacobi(n,1000,tol,A,V,time,count);
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
