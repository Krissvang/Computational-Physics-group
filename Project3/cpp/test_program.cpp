#include <iostream>
#include <cmath>
#include "lib.h"

using namespace std;

int main(int nchar, char* argv[])
{
  long idum = -1;
  for (int i = 0; i < 5; i++)
  {
     cout << 2*acos(-1.)*ran0(&idum) <<endl;
    /* code */
  }
  


  return 0;
}