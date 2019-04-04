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
  // Message
  double x_min, x_max, y_min, y_max, t , tf , dt, test,hy,hx;
  int Nx,Ny,n;
  test = 10;
  x_min = -0.005;
  x_max = 0.005;
  y_min = -0.0025;
  y_max = 0.0025;
  Nx = 250;
  Ny = 125;
  t = 0.;
  tf = 100;
  dt = 1*pow(10,-1);

  hy = (y_max - y_min) / Ny;
  hx = (x_max - x_min) / Nx;

  n= (tf-t)/dt;

  Matrices* sys(0);

  cout << 1 << endl;
  sys = new Matrices(x_min,x_max,y_min,y_max,Nx,Ny);
  cout << "1" << endl;
  sys -> Rho(0);
  cout << "Rho ok" << endl;
  sys -> Xi();
  cout << "Xi ok" << endl;
  sys -> Lambda();
  cout << "Lambda ok" << endl;
  sys -> A();
  cout << "A ok" << endl;
  sys -> AStar();
  cout << "Astar ok" << endl;
  sys -> L1234();
  cout << "L1234 ok" << endl;
  sys -> R();
  cout << "R ok" << endl;
  sys -> M();
  cout << "M ok" << endl;
  sys -> Flux(t);
  cout << "F ok" << endl;


  // cout << "Rho" << endl;
  // cout << sys->GetRho()<<endl;
  // cout << "Xi" << endl;
  // cout << sys->GetXi() << endl;
  // cout << "Lambda" << endl;
  // cout << sys->GetLambda() << endl;
  // cout << "A" << endl;
  // cout << sys->GetA() << endl;
  // cout << "L1" << endl;
  // cout << sys->GetL1() << endl;
  // cout << "L2" << endl;
  // cout << sys->GetL2() << endl;
  //cout << "L3" << endl;
  //cout << sys->GetL3() << endl;
  //cout << "L4" << endl;
  //cout << sys->GetL4() << endl;  _hy = (_y_max - _y_min) / _Nc;


  //cout << "R" << endl;
  //cout << sys->GetR() << endl;

  //cout << "M" <<endl;
//  cout << sys->GetM() << endl;

  cout << "Flux" << endl;
  cout << sys->GetFlux();
  int i;
  i=0;

  ofstream sol;
  sol.open("solution.dat",ios::out);

  while(t<tf) //faire condition while sur dt ensuite
{
//cout << "t="<<t << endl;
sys->RhoStar(t);
sys->Flux(t);
//cout << sys -> GetFlux() << endl;
//cout << "Newton" <<endl;
sys->iteration();
//sys->SaveSolPara("Solution");
sol  << t << " " << sys->GetTvec().coeffRef(249 * Ny+75) << " " << 249 * hx << " " << 75 * hy << " " << sys->GetRho().coeffRef(0,75) << endl;
cout << t << endl;

//cout << "Rho" <<endl;
sys->Rho(t);
t+=dt;
i+=1;

}
sol.close();


  return 0;
}
