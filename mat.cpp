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
//#include <Eigen/SparseCholesky>

using namespace Eigen;

using namespace std;

Matrices::Matrices( const double x_min, const double x_max, const double y_min, const double y_max, const int Nx , const int Ny):
_Nl(Nx),_Nc(Ny),_x_min(x_min), _x_max(x_max), _y_min(y_min), _y_max(y_max)
{
  // Paramètres
  _Ar = 1000;
  _M.resize(_Nl*_Nc,_Nl*_Nc);
  _Lm = 300000;
  _hy = (_y_max - _y_min) / _Nc;
  _hx = (_x_max - _x_min) / _Nl;
  _RhoV = 1500;
  _RhoP = 1000;
  _Sm.resize(_Nl*_Nc);
  _Rho.resize(_Nl,_Nc);
  _RhoStar.resize(_Nl,_Nc);
  _Xi.resize(_Nl,_Nc);
  _XiStar.resize(_Nl,_Nc);
  _Lambda.resize(_Nl,_Nc);
  _A.resize(_Nl,_Nc);
  _AStar.resize(_Nl,_Nc);
  _T.resize(_Nl,_Nc);

  _Tvect.resize(_Nl*_Nc);
  _R.resize(_Nl,_Nc);
  _LambV = 1;
  _LambP = 1;
  _Cpp = 1500;
  _Cpv = 1000;
  _T0 = 293;
  _TA = 6000;
  _dt = pow(10,-2) ;
  _t = 0;
  _T.setZero();
  _Tvect.setZero();
  _Rho.setZero();
  _Lambda.setZero();
  _L1.resize(_Nl,_Nc); _L2.resize(_Nl,_Nc); _L3.resize(_Nl,_Nc); _L4.resize(_Nl,_Nc);
  _FS.resize(_Nc); _FN.resize(_Nc); _FO.resize(_Nl); _FE.resize(_Nl);
  _f.resize(_Nc*_Nl);
  //_P1.resize(_Nc,_Nc); _P2.resize(_Nc,_Nc); _P3.resize(_Nc,_Nc); _P4.resize(_Nc,_Nc);

  for (int i=0 ; i<_Nl; i++)
  {
    //cout << "i = " << i << endl;
    for (int j=0 ; j<_Nc ; j++)
    {
      _T.coeffRef(i,j) = _T0;
      _Tvect.coeffRef(j+_Nc*i) = _T.coeffRef(i,j);
      _Rho.coeffRef(i,j) = _RhoV;
      _Lambda.coeffRef(i,j) = 1;
      _A.coeffRef(i,j) = _Cpv;
      _AStar.coeffRef(i,j) = _Cpv;
    }
  }
}

void Matrices::Flux(double t)
{
  _FS.setZero(); _FN.setZero(); _FO.setZero(); _FE.setZero(); _f.setZero();

  for (int i=0; i<_FS.size(); i++)
  {
    _FS(i) = 0;
    if (t<=50)
    {
      _FN(i) = 10000*t;
    }
    else
    {
      _FN(i) = 500000 - 9000 * (t-50);
    }
  }

  for (int i=0 ; i<_FO.size(); i++)
  {
    _FO(i) = 0;
    _FE(i) = 0;
  }

  for (int j=0; j<_Nc; j++)
  {
    //cout << i <<endl;
    _f(j) += _dt / ( _A(1,j)* _hx) * _FS(j);
    _f(_Nl*_Nc - 1 -j) += _dt / (_A(_Nl-1,j) * _hx) * _FN(_Nc - 1 - j);
  }

  for (int i=0; i<_Nl; i++)
  {
    _f(i*_Nc) +=_dt / (_A(i,1)*_hy) * _FO(i);
    _f((i+1)*_Nc -1) += _dt / (_A(i,_Nc-1)*_hy) * _FE(i);
  }
}


void Matrices::Rho(double t)
{
  double C;


  for (int i=0; i<_Nl; i++)
  {
    for (int j=0; j<_Nc ; j++)
    {

      C = _RhoV * _Ar * exp(-_TA/_Tvect.coeffRef(i*_Nc+j)) / (_RhoV-_RhoP);
      _Rho.coeffRef(i,j) = (_RhoV-_RhoP) * exp(-C*t) + _RhoP;

    }
  }
}

void Matrices::RhoStar(double t)
{
  double C;
  for (int i=0; i<_Nl; i++)
  {
    for (int j=0; j<_Nc ; j++)
    {
      C = _RhoV * _Ar * exp(-_TA/_Tvect.coeffRef(i*_Nc+j)) / (_RhoV-_RhoP);
      _RhoStar.coeffRef(i,j) = (_RhoV-_RhoP) * exp(-C*(t)) + _RhoP;

    }
  }
}

void Matrices::Xi()
{
  for (int i=0; i<_Nl; i++)
  {
    for (int j=0 ; j<_Nc; j++)
    {
      _Xi.coeffRef(i,j) = (_RhoV - _Rho.coeffRef(i,j)) / (_RhoV-_RhoP);
    }
  }
}

void Matrices::XiStar()
{
  for (int i=0; i<_Nl; i++)
  {
    for (int j=0 ; j<_Nc; j++)
    {
      _XiStar.coeffRef(i,j) = (_RhoV - _RhoStar.coeffRef(i,j)) / (_RhoV-_RhoP);
    }
  }
}

void Matrices::Lambda()
{
  for (int i=0; i<_Nl; i++)
  {
    for (int j=0 ; j<_Nc ; j++)
    {
      _Lambda.coeffRef(i,j) = (1- _Xi.coeffRef(i,j)) * _LambV + _Xi.coeffRef(i,j) * _LambP ;
    }
  }
}

void Matrices::A()
{
  for (int i=0 ; i<_Nl; i++)
  {
    for (int j=0 ; j<_Nc ; j++)
    {
      _A.coeffRef(i,j) = (1 - _Xi.coeffRef(i,j)) * _Cpv * _RhoV + _Xi.coeffRef(i,j) * _Cpp * _RhoP;
    }
  }
}
void Matrices::AStar()
{
  for (int i=0 ; i<_Nl; i++)
  {
    for (int j=0 ; j<_Nc ; j++)
    {
      _AStar.coeffRef(i,j) = (1 - _XiStar.coeffRef(i,j)) * _Cpv * _RhoV + _XiStar.coeffRef(i,j) * _Cpp * _RhoP;
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
      if (_A(i,j)!=0)
      {
        _R(i,j) = _AStar(i,j) / _A(i,j);
      }
      else
      {
        cout << "erreur dans R, A=0" << endl;
      }
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
      triplets.push_back({i*_Nc+j,i*_Nc+j,_R(i,j)+(_L1(i,j)+_L2(i,j)+_L3(i,j)+_L4(i,j))});
      if(j+_Nc*i<_Nl*_Nc-1)
      {
        triplets.push_back({i*_Nc+j,i*_Nc+j+1,-_L2(i,j)});
      }
      if(i<_Nl-1)
      {
        triplets.push_back({i*_Nc+j,(i+1)*_Nc+j,-_L4(i,j)});
      }
      if (j+_Nc*i>0)
      {
        triplets.push_back({i*_Nc+j,i*_Nc+j-1,-_L1(i,j)});
      }
      if (i>0)
      {
        triplets.push_back({i*_Nc+j,i*_Nc+j-_Nc,-_L3(i,j)});
      }
    }
  }
  _M.setFromTriplets(triplets.begin(),triplets.end());
}

void Matrices::Sm()
{
  for (int i=0; i<_Nl; i++)
  {
    for (int j=0; j<_Nc; j++)
    {
      _Sm(i*_Nc+j) = _T0 * (1-_R(i,j)) + (_Lm/_A(i,j))* (_Rho(i,j)-_RhoStar(i,j));
    }
  }
}

void Matrices::iteration()
{
  VectorXd sol1;
  VectorXd zero;
  zero.resize(_Tvect.size());
  for (int i=0 ;  i<_Tvect.size(); i++)
  {
    zero(i)=0.;
  }
  //cout << "début BiCGSTAB" << endl;
  BiCGSTAB <SparseMatrix<double>> solver;

  //cout<< "Mise à jour données" << endl;
  XiStar(); Xi(); AStar(); A();
  Lambda(); L1234(); R(); Sm(); M();// Flux();

  //cout << "Sm" <<_Sm.norm() << endl;
  solver.compute(_M);
  //cout << "compute ok" << endl;

  sol1= solver.solve(_Tvect+_f-_Sm);
  //cout << "itérations = " << solver.iterations() << endl;
  //cout << "erreur estimée = " << solver.error() << endl;


  _Tvect=sol1;



}

void Matrices::SaveSol(ofstream sol, int i, int j, double t)
{
      sol  << t << " " << _Tvect(i * _Nc+j) << " " << i * _hx << " " << j * _hy << " " << _Rho(i,j) << endl;
}

void Matrices::SaveSolPara(string nf)
{
  ofstream solution;
  solution.open(nf, ios::out);
  solution.precision(7);
  solution << "# vtk DataFile Version 3.0" << endl;
  solution << "sol" << endl;
  solution << "ASCII" << endl;
  solution << "DATASET STRUCTURED_POINTS" << endl;
  solution << "DIMENSIONS " << _Nl << " " << _Nc << " " << 1 << endl;
  solution << "ORIGIN " << _x_min << " " << _y_min << " " << 0 << endl;
  solution << "SPACING " << _hx << " " << _hy << " " << 1 << endl;;
  solution << "POINT_DATA " << _Nc*_Nl << endl;
  solution << "SCALARS sol float" << endl;
  solution << "LOOKUP_TABLE default" << endl;
  for(int i=0; i<_Nl; i++)
  {
    for(int j=0; j<_Nc; j++)
    {
      solution << float(_Tvect(i*_Nc+j)) << " ";
    }
    solution << endl;
  }
  solution.close();
}
