#ifndef JACOBI_H
#define	JACOBI_H

#include <iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include <armadillo>
#include<time.h>
#include<vector>

using namespace std;
using namespace arma;

void print_vals(mat,mat,int,double);
void initialize(int,double,mat&,vec&,mat&,int,double);
int jacobi(int, int, double, mat& , mat& , double& , int& );
void find_max(mat,int&,int&,double&,int);
vector<double> get_eigenvals(mat,int);
mat get_eigenvecs(mat,mat,int);
void initialize_classic(int, double, mat& , vec& , mat& ,int ,double );
void initialize_schrodinger(int, double, mat&,vec&, mat&,int,double);
void initialize_beam(int , double , mat& , vec& , mat&, int ,double );
void write_eigenpairs(vec, mat, vec, double, double, string, int,int);
void output_count_beam(int , string, double, double, int);

#endif /* JACOBI_H */
