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
#include "pyro.h"
//#include <Eigen/SparseCholesky>

using namespace Eigen;

using namespace std;

Matrices::Matrices( const double x_min, const double x_max, const double y_min, const double y_max):
_Nx(Nx),_Ny(Ny),_x_min(x_min), _x_max(x_max), _y_min(y_min), _y_max(y_max)
{
  // Paramètres
  _M.resize(_Nx*_Ny,_Nx*_Ny);
  _Lm = 10;
  _RhoV = 1500;
  _RhoP = 2;
  _Rho.resize(_Nx,_Ny);
  _RhoStar.resize(_Nx,_Ny);
  _Xi.resize(_Nx,_Ny);
  _XiStar.resize(_Nx,_Ny);
  _Lambda.resize(_Nx,_Ny);
  _A.resize(_Nx,_Ny);
  _AStar.resize(_Nx,_Ny);
  _T.resize(_Nx,_Ny);
  _LambV = 3;
  _LambP = 4;
  _Cpp = 1000;
  _Cpv = 6;
  _T0 = 293;
  _dt = 20 ;
  _t = 0;
  _L1.resize(_Nx,_Ny); _L2.resize(_Nx,_Ny); _L3.resize(_Nx,_Ny); _L4.resize(_Nx,_Ny);
  _P1.resize(_Nx,_Nx); _P2.resize(_Nx,_Nx); _P3.resize(_Nx,_Nx); _P4.resize(_Nx,_Nx);

  for (int i=0 ; i<_Nx; i++)
  {
    cout << "i = " << i << endl;
    for (int j=0 ; j<_Ny ; j++)
    {
      cout << "j = " << j << endl;
      _T(i,j) = _T0;
      _Sol(j+_Ny*i) = _T(i,j);
      _Rho(i,j) = _RhoV;
      _Lambda(i,j) = 1;
      _A(i,j) = Cpv;
      _AStar(i,j) = Cpv;
    }
  }
}

void Matrices::Xi()
{
  for (int i=0; i<_Nx; i++)
  {
    for (int j=0 ; j<_Ny; j++)
    {
      _Xi(i,j) = (_RhoV - _Rho(i,j)) / (_RhoV-_RhoP);
    }
  }
}

void Matrices::Lambda()
{
  for (int i=0; i<_Nx; i++)
  {
    for (int j=0 ; j<_Ny ; j++)
    {
      _Lambda(i,j) = (1- _Xi(i,j)) * _LambV + _Xi(i,j) * _LambP ;
    }
  }
}

void Matrices::A()
{
  for (int i=0 ; i<_Nx; i++)
  {
    for (int j=0 ; j<_Ny ; j++)
    {
      _A(i,j) = (1 - _Xi(i,j)) * _Cpv * _RhoV + _Xi(i,j) * _Cpp * _RhoP;
    }
  }
}

void Matrices::L1234()
{
  _L1.setZero();_L2.setZero();_L3.setZero();_L4.setZero();
  for (int i=0; i<_Ny; i++)
  {
    for (int j=1 ; j<_Nx; j++)
    {
      _L1(i,j) = _dt * (_Lambda(i,j)+_Lambda(i,j-1))/(2*_hy*_hy*_A(i,j));
      _L2(i,j-1) = _dt * (_Lambda(i,j-1)+_Lambda(i,j))/(2*_hy*_hy*_A(i,j-1));
    }
  }
  for (int i=1; i<_Ny; i++)
  {
    for (int j=0; j<_Nx ; j++)
    {
      _L3(i,j) = _dt * (_Lambda(i,j)+_Lambda(i-1,j))/(2*_hx*_hx*_A(i,j));
      _L4(i-1,j) = _dt * (_Lambda(i-1,j)+_Lambda(i,j))/(2*_hx*_hx*_A(i-1,j));
    }
  }
}

void Matrices::R()
{
  for (int i=0; i<_Ny; i++)
  {
    for (int j=0; j<_Nx; j++)
    {
      R(i,j) = AStar(i,j) / A(i,j);
    }
  }
}

void Matrices::PP14()
{
  //Définition de P4, y=0, x varie
  vector<Triplet<double>> triplets4;
  triplets4.push_back({0,0,_R(0,0)-(_L2(0,0)+_L4(0,0))});
  triplets4.push_back({1,0,_L2(0,1)});
  //Définition de P1, y= y_max, x varie
  vector<Triplet<double>> triplets1;
  triplets3.push_back({0, 0,_R(0,_Ny-1)-(2*_L1(0,_Ny-1)+ _L2(0,_Ny-1))});
  triplets3.push_back({1, 0,_L1(1,_Ny-1)});

  for(int i=1 ; i<_Nx-1; i++)
  {
    triplets4.push_back({i,i,_R(i,0)-(_L1(i,0)+_L2(i,0)+_L1(i,0))});
    triplets4.push_back({i-1,i,_L2(i-1,0)});
    triplets4.push_back({i+1,i,_L2(i+1,0)});

    triplets1.push_back({i,i,_R(i,_Ny-1)-(_L1(i,_Ny-1)+_L2(i,_Ny-1)+_L1(i,_Ny-1))});
    triplets1.push_back({i+1,i,_L1(i+1,_Ny-1)}});
    triplets1.push_back({i-1,i,_L1(i-1,_Ny-1)}});
  }

  triplets1.push_back({_Nx-1,_Nx-1,_R(_Nx-1,_Ny-1)-(2*_L1(_Nx-1,_Ny-1)+ _L3(_Nx-1,_Ny-1))});
  triplets1.push_back({_Nx-2,_Nx-1, _L1(_Nx-2,_Ny-1)});
  triplets4.push_back({_Nx-1,_Nx-1, _R(_Nx-1,0)-(2*_L2(_Nx-1,0)+ _L3(_Nx-1,0)) });
  triplets4.push_back({_Nx-2,_Nx-1, _L2(_Nx-2,0)});

  P1.setFromTriplets(triplets1.begin(),triplets1.end());
  P4.setFromTriplets(triplets4.begin(),triplets4.end());


}

void Matrices::M()
{

}
