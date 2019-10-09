//   This is a simple program which tests Gaussian quadrature using Legendre and Laguerre polynomials
//   It integrates the simple function x* exp(-x) for the interval
//   x \in [0,infty). The exact result is 1. For Legendre based quadrature a
//   tangent mapping is also used.

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#define EPS 3.0e-14
#define MAXIT 10
#define   ZERO       1.0E-10
using namespace std;

//     Here we define various functions called by the main program

double int_function(double x1, double y1, double z1, double x2, double y2, double z2);
double one_electron_function(double x);
void gauleg(double x1, double x2, double x[], double w[], int n);


//   Main function begins here
int main()
{
    ofstream ofile;
    string outfile;
    int n;
    double a, b;
    int selection;
    
    cout << "Insert 1 if you want to evaluate the one electron function, insert 0 if you want to solve the integral" << endl;
    cin >> selection;
    
    if(selection==1){
        cout << "Read in the number of points" << endl;
        cin >> n;
        cout << "Read the function limits" << endl;
        cin >> a >> b;
        cout << "Write the output file name" <<endl;
        cin >> outfile;
        double *x = new double [n];
        double *f = new double [n];
        double h=(b-a)/n;
    
        for (int i=0; i<n; i++){
            x[i]=a+i*h;
        }

        for(int i=0; i<n; i++){
            f[i]=one_electron_function(x[i]);
        }
    
        ofile.open(outfile);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << "number of points:" << endl;
        ofile << n << endl;
        ofile << "      points:    One electron function(points):" << endl;
        for(int i=0; i<n; i++){
            ofile << setw(19) << setprecision(8) << x[i];
            ofile << setw(19) << setprecision(8) << f[i] <<endl;
        }
        
        delete [] x;
        delete [] f;
    }
    
    if(selection==0){
        int N;
        cout << "Read in the number of integration points" << endl;
        cin >> N;
        cout << "Read in integration limits" << endl;
        cin >> a >> b;
        
        //reserve space in memory for vectors containing the mesh points
        //weights and function values for the use of the gauss-legendre
        //method
        double *x = new double [N];
        double *w = new double [N];
        //set up the mesh points and weights
        gauleg(a,b,x,w, N);

        //evaluate the integral with the Gauss-Legendre method
        //Note that we initialize the sum
        double int_gauss = 0.;
        //six-double loops
        for (int i=0;i<N;i++){
            for (int j = 0;j<N;j++){
                for (int k = 0;k<N;k++){
                    for (int l = 0;l<N;l++){
                        for (int m = 0;m<N;m++){
                            for (int n = 0;n<N;n++){
                                int_gauss+=w[i]*w[j]*w[k]*w[l]*w[m]*w[n]
                                    *int_function(x[i],x[j],x[k],x[l],x[m],x[n]);
                    }}}}}
            }
        cout << int_gauss << endl;
    }
    
    return 0;
}  // end of main program

//  this function defines the function to integrate
double int_function(double x1, double y1, double z1, double x2, double y2, double z2)
{
   double alpha = 2.;
// evaluate the different terms of the exponential
    double exp1=-2*alpha*sqrt(x1*x1+y1*y1+z1*z1);
    double exp2=-2*alpha*sqrt(x2*x2+y2*y2+z2*z2);
    double deno=sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2));
    if(deno==0){deno=1e-10;}
    return exp(exp1+exp2)/deno;
} // end of function to evaluate



//  this function defines the wave function of one electron
double one_electron_function(double x)
{
    double value = exp(-2*x);
    return value;
} // end of function to evaluate


       /*
       ** The function
       **              gauleg()
       ** takes the lower and upper limits of integration x1, x2, calculates
       ** and return the abcissas in x[0,...,n - 1] and the weights in w[0,...,n - 1]
       ** of length n of the Gauss--Legendre n--point quadrature formulae.
       */

void gauleg(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
       ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
     p2 =0.0;

          /*
       ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

     for(j = 1; j <= n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
     }

       /*
       ** p1 is now the desired Legrendre polynomial. Next compute
           ** ppp its derivative by standard relation involving also p2,
           ** polynomial of one lower order.
           */
 
     pp = n * (z * p1 - p2)/(z * z - 1.0);
     z1 = z;
     z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
      ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
} // End_ function gauleg()
