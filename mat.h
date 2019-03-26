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

class Matrices{
private:
  const double _x_min;
  const double _x_max;
  const double _y_min;
  const double _y_max;
  const int _Nx;
  const int _Ny;
  const double _RhoV,_RhoP,_LambV,_LambP;
  const double _Cpp , Cpv, _Lm;
  double _L1 , _L2 , _L3 , _L4;
  double _t, _dt;
  double _hx,_hy;
  double _T0;

  Eigen::SparseMatrix<double> _M;
  Eigen::VectorXd _x;
  Eigen::VectorXd _y;
  Eigen::VectorXd _f;
  Eigen::VectorXd _Sol;
  Eigen::MatrixXd _L1 , _L2 , _L3 , _L4, , _R;
  Eigen::MatrixXd _T , _A, _AStar , _Lambda , _RhoStar , _Xi , _XiStar;
  Eigen::SparseMatrix<double> _P, _P1 , _P4 ;
  Eigen::VectorXd _AVec , _LambdaVec , _RhoVec , _XiVec ;
  Eigen::VectorXd _FS , _FN , _FE , _FO;

public:
Matrices(const double x_min, const double x_max, const double y_min, const double y_max,const int Nx, const int Ny);

void Xi();
void Lambda();
void A();
void M();
void L1234();
void R();
}

};
