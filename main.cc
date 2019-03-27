#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include "Dense"
#include "Sparse"
#include <chrono>
#include "mat.h"

using namespace Eigen;

using namespace std;

int main() {
  // Je fais un test
  double x_min, x_max, y_min, y_max, t , tf , dt, test;
  int Nx,Ny,n;
  test = 10;
  x_min = -0.0025;
  x_max = 0.0025;
  y_min = -0.005;
  y_max = 0.005;
  Nx = 125;
  Ny = 250;
  t = 0.;
  tf = 100;
  dt = pow(10,-1);

  n= (tf-t)/dt;

  Matrices* sys(0);

  cout << 1 << endl;
  sys = new Matrices(x_min,x_max,y_min,y_max,Nx,Ny);

  sys -> Rho(0);
  cout << sys->GetRho() << endl;

  return 0;
}
