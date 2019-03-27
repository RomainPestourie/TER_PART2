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
  const int _Nl;
  const int _Nc;
  double _RhoV,_RhoP,_LambV,_LambP;
  double _Cpp , _Cpv, _Lm, _Ar;
  double _t, _dt;
  double _hx,_hy;
  double _T0,_TA;

  Eigen::SparseMatrix<double> _M,_T;
  Eigen::VectorXd _x;
  Eigen::VectorXd _y;
  Eigen::VectorXd _f;
  Eigen::VectorXd _Sol;
  Eigen::MatrixXd _L1 , _L2 , _L3 , _L4 , _R;
  Eigen::MatrixXd _Rho , _A, _AStar , _Lambda , _RhoStar , _Xi , _XiStar;
  Eigen::VectorXd _Tvec ;
  Eigen::VectorXd _FS , _FN , _FE , _FO;

public:
Matrices(const double x_min, const double x_max, const double y_min, const double y_max,const int Nx, const int Ny);

void Rho(double t);
void Xi();
void Lambda();
void A();
void M();
void L1234();
void R();
Eigen::VectorXd GetT() {return _T;};
Eigen::SparseMatrix<double> GetM() {return _M;};
Eigen::MatrixXd GetRho() {return _Rho;};
Eigen::MatrixXd GetA() {return _A;};
Eigen::MatrixXd GetXi() {return _Xi;};
Eigen::MatrixXd GetLambda() {return _Lambda;};
};
