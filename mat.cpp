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
_Nx(Nl),_Ny(Nc),_x_min(x_min), _x_max(x_max), _y_min(y_min), _y_max(y_max)
{
  // Paramètres
  _M.resize(_Nl*_Nc,_Nl*_Nc);
  _Lm = 10;
  _RhoV = 1500;
  _RhoP = 2;
  _Rho.resize(_Nl,_Nc);
  _RhoStar.resize(_Nl,_Nc);
  _Xi.resize(_Nl,_Nc);
  _XiStar.resize(_Nl,_Nc);
  _Lambda.resize(_Nl,_Nc);
  _A.resize(_Nl,_Nc);
  _AStar.resize(_Nl,_Nc);
  _T.resize(_Nl,_Nc);
  _LambV = 3;
  _LambP = 4;
  _Cpp = 1000;
  _Cpv = 6;
  _T0 = 293;
  _dt = 20 ;
  _t = 0;
  _L1.resize(_Nl,_Nc); _L2.resize(_Nl,_Nc); _L3.resize(_Nl,_Nc); _L4.resize(_Nl,_Nc);
  //_P1.resize(_Nc,_Nc); _P2.resize(_Nc,_Nc); _P3.resize(_Nc,_Nc); _P4.resize(_Nc,_Nc);
  _P.resize(_Nc*_Nl,_Nc*_Nl);

  for (int i=0 ; i<_Nl; i++)
  {
    cout << "i = " << i << endl;
    for (int j=0 ; j<_Nc ; j++)
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
  for (int i=0; i<_Nl; i++)
  {
    for (int j=0 ; j<_Nc; j++)
    {
      _Xi(i,j) = (_RhoV - _Rho(i,j)) / (_RhoV-_RhoP);
    }
  }
}

void Matrices::Lambda()
{
  for (int i=0; i<_Nl; i++)
  {
    for (int j=0 ; j<_Nc ; j++)
    {
      _Lambda(i,j) = (1- _Xi(i,j)) * _LambV + _Xi(i,j) * _LambP ;
    }
  }
}

void Matrices::A()
{
  for (int i=0 ; i<_Nl; i++)
  {
    for (int j=0 ; j<_Nc ; j++)
    {triplets.push_back({0,0,_R(0,0)-(L1(0,0)+_L2(0,0)+L3(0,0)+_L4(0,0))});
  triplets.push_back({1,0,_L2(0,1)});
      _A(i,j) = (1 - _Xi(i,j)) * _Cpv * _RhoV + _Xi(i,j) * _Cpp * _RhoP;
    }
  }
}

void Matrices::L1234()
{
  _L1.setZero();_L2.setZero();_L3.setZero();_L4.setZero();
  for (int i=0; i<_Nl; i++)
  {
    for (int j=1 ; j<_Nc; j++)
    {
      _L1(i,j) = _dt * (_Lambda(i,j)+_Lambda(i,j-1))/(2*_hy*_hy*_A(i,j));
      _L2(i,j-1) = _dt * (_Lambda(i,j-1)+_Lambda(i,j))/(2*_hy*_hy*_A(i,j-1));
    }
  }
  for (int i=1; i<_Nl; i++)
  {
    for (int j=0; j<_Nc ; j++)
    {
      _L3(i,j) = _dt * (_Lambda(i,j)+_Lambda(i-1,j))/(2*_hx*_hx*_A(i,j));
      _L4(i-1,j) = _dt * (_Lambda(i-1,j)+_Lambda(i,j))/(2*_hx*_hx*_A(i-1,j));
    }
  }
}

void Matrices::R()
{
  for (int i=0; i<_Nl; i++)
  {
    for (int j=0; j<_Nc; j++)
    {
      R(i,j) = AStar(i,j) / A(i,j);
    }
  }
}

void Matrices::M()
{
  vector<Triplet<double>> triplets;

  for (int i=0 ; i<_Nl; i++)
  {
    for (int j=0 ; j<_Nc;j++)
    {
      triplets.push_back({i*_Nc+j,i*_Nc+j,R(i,j)-(_L1(i,j)+_L2(i,j)+_L3(i,j)+_L4(i,j))})
      if (i*_Nc+j>1)
      {
        triplets.push_back({i*_Nc+j-1,i*_Nc+j,_L1(i,j)})
      }
      else if(i*_Nc+j+1<(_Nl-1)*(_Nc-1))
      {
        triplets.push_back({i*_Nc+j+1,i*_Nc+j,_L2(i,j)})
      }
    }
  }





  }
}

void Matrices::M()
{

}
