#include <iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<chrono>
#include<algorithm>
#include<vector>

#include "jacobi.h"

using namespace std::chrono;

// performs jacobi algorithm
// to find eigenvalues/vectors
int jacobi(int n, int maxcount, double conv, mat& a, mat& v, double& time, int& count) {
    cout.precision(5);
    double aip=0, aiq=0, vip=0, viq=0;
    double tau=0, t=0, s=0, c=0;//tan(theta), sin(theta), cos(theta)
    count=1;                    //count of iterations
    int count_old=count-10;     //keep track of every 10th iteration
    int p=n-1, q=n-2;           //off diag all same value to start
                                //pick last as first maximum
    time_point<high_resolution_clock> start, end;

    if(n<=10){
        cout<<"Before diagonalization"<<endl;
        print_vals(a,v,n,conv);
        cout<<endl;
    }

    double app=a(p,p);
    double aqq=a(q,q);
    double apq=a(p,q);

    start = high_resolution_clock::now();

    while(abs(apq)>conv){
        if(count>1){
            apq=0;
            find_max(a,p,q,apq,n);
        }

        //calculate sin(theta) and cos(theta)
        aqq=a(q,q);
        app=a(p,p);
        tau=(aqq-app)/(2*apq);
        if(tau>0)
            t=1./(tau+sqrt(1+tau*tau));
        else
            t=-1./(-tau+sqrt(1+tau*tau));
        c=1/sqrt(1+t*t);
        s=c*t;

        //calculate new matrix elements and vectors
        for(int i=0;i<n;i++){
            if(i!=p && i!=q){
                aip=a(i,p);
                aiq=a(i,q);
                a(i,p)=aip*c-aiq*s;
                a(p,i)=aip*c-aiq*s;
                a(i,q)=aiq*c+aip*s;
                a(q,i)=aiq*c+aip*s;
            }
            vip=v(i,p);
            viq=v(i,q);
            v(i,p)=c*vip-s*viq;
            v(i,q)=c*viq+s*vip;
        }
        a(p,p)=app*c*c-2*apq*c*s+aqq*s*s;
        a(q,q)=app*s*s+2*apq*c*s+aqq*c*c;
        a(p,q)=0;
        a(q,p)=0;

        count++;
    }

    end=high_resolution_clock::now();

    if(n<=10){
        cout<<"After diagonalization"<<endl;
        print_vals(a,v,n,conv);
        cout<<endl;
    }

    cout<<"Diagonalization took "<<count<<" iterations"<<endl;
    duration<double> elapsed = end-start;
    time = elapsed.count();
    cout<<scientific<<"CPU time (sec) : "<<time<<endl;

    return 0;
}

//get first three eigenvectors
mat get_eigenvecs(mat a, mat v, int n){
    vector<double>eigenvals=get_eigenvals(a,n);
    mat vecs(3,n);
    for(int i=0;i<3;i++){
        for(int j=0;j<n;j++){
            if(a(j,j)==eigenvals[i]){
                for(int k=0;k<n;k++){
                      vecs(i,k)=v(k,j);
                }
             }
         }
    }
    return vecs;
}

//get eigenvalues in order
vector<double> get_eigenvals(mat a,int n){
    vector<double>eigen;
    for(int i=0;i<n;i++){
        eigen.push_back(a(i,i));
    }
    sort (eigen.begin(), eigen.begin()+n);
    return eigen;
}

//initialize matrix/vectors
void initialize_schrodinger(int n, double h, mat& a, vec& r, mat& v,int interact,double wr){
    //initialize x values
    r(0)=h;
    for (int i=1; i<n ;i++){
        r(i)=r(i-1)+h;
    }
    double h2inv = 1/(h*h);
    //initialize matrix and vector
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j && interact==0){
                a(i,j)=2*h2inv+r(i)*r(i);
                v(i,j)=1;
            }
            else if (i==j && interact==1){
                a(i,j)=2*h2inv+wr*wr*r(i)*r(i)+1/r(i);
                v(i,j)=1;
            }
            else if (i==j+1 or i==j-1){
                a(i,j)=-1*h2inv;
            }
            else{
                a(i,j)=0;
                v(i,j)=0;
            }
        }
    }
}
//initialize matrix/vectors
void initialize_beam(int n, double h, mat& a, vec& r, mat& v,int interact,double wr){
    //initialize x values
    r(0)=h;
    for (int i=1; i<n ;i++){
        r(i)=r(i-1)+h;
    }
    double h2inv = 1/(h*h);
    //initialize matrix and vector
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j && interact==0){
                a(i,j)=2*h2inv;
                v(i,j)=1;
            }
            else if (i==j && interact==1){
                a(i,j)=2*h2inv;
                v(i,j)=1;
            }
            else if (i==j+1 or i==j-1){
                a(i,j)=-1*h2inv;
            }
            else{
                a(i,j)=0;
                v(i,j)=0;
            }
        }
    }
}

//initialize matrix/vectors
void initialize_classic(int n, double h, mat& a, vec& r, mat& v,int interact,double wr){
    //initialize x values
    r(0)=h;
    for (int i=1; i<n ;i++){
        r(i)=r(i-1)+h;
    }

    //initialize matrix and vector
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j && interact==0){
                a(i,j)=2;
                v(i,j)=1;
            }
            else if (i==j && interact==1){
                a(i,j)=2;
                v(i,j)=1;
            }
            else if (i==j+1 or i==j-1){
                a(i,j)=-1;
            }
            else{
                a(i,j)=0;
                v(i,j)=0;
            }
        }
    }
}


//find maximum non-diag matrix elements
void find_max(mat a,int& p,int& q,double& apq,int n){
    for (int i=0;i<n;i++){
         for (int j=0;j<n;j++){
            if(i!=j && abs(a(i,j))>=abs(apq)){
                apq=a(i,j);
                p=i;
                q=j;
            }
         }
    }
}

//print matrix and eigenvectors
void print_vals(mat A, mat v,int n,double conv){
    cout<<"A: ";
    for (int i=0;i<n;i++){
        if(i>0){
            cout<<"   ";
         }
        for (int j=0;j<n;j++){
            if(abs(A(i,j))>conv)
                cout<<fixed<<A(i,j)<<" ";
            else cout<<"0.000 ";
        }
        cout<<endl;
    }
    for (int i=0;i<n;i++){
        cout<<"v"<<i<<": ";
        for (int j=0;j<n;j++){
            if(abs(v(j,i))>conv)
                cout<<fixed<<v(j,i)<<" ";
            else cout<<"0.000 ";
        }
        cout<<endl;
    }
}

void write_eigenpairs(vec jacobi_e_vals, mat e_vecs, vec armadillo_eigenval, double t_Jacobi, double t_arma, string outfile,int count, int n){
    outfile="output/"+outfile;
    ofstream ofile;
    ofile.open(outfile);
    ofile << setiosflags(ios::showpoint);
    ofile << "t_Jacobi: " << scientific << t_Jacobi << " (" << count
    << " iterations) t_armadillo: " << scientific << t_arma << endl;
    ofile << "Jacobi eigenvalues:    Armadillo eigenvalues:                v0:"
    "                v1:                v2:" << endl;
    for(int i = 0; i < n; i++){
      ofile << setw(19) << setprecision(8) << jacobi_e_vals[i];
      ofile << setw(26) << setprecision(8) << armadillo_eigenval[i];
      ofile << setw(19) << setprecision(8) << e_vecs(0,i);
      ofile << setw(19) << setprecision(8) << e_vecs(1,i);
      ofile << setw(19) << setprecision(8) << e_vecs(2,i) << endl;
    }
    ofile.close();
}

void output_count_beam(int n, string filename, double h, double wr, int interact ){
    ofstream ofile;
    filename="output/"+filename;
    ofile.open(filename);
    ofile << setiosflags(ios::showpoint);
    ofile << "  Count:       CPU time (s):"<< endl;
    double time=0;
    for (int i = 2; i < n; i++)
    {
        int count=0;
        mat A_temp(i,i);
        mat V_temp(i,i);
        vec r_temp(i);
        initialize_beam(i, h, A_temp, r_temp, V_temp, interact, wr);
        jacobi(i,500,0.0000000001,A_temp,V_temp,time,count);
        ofile << setw(8) << count;
        ofile << setw(20) << setprecision(8) << time << endl;
    }
    ofile.close();
};
