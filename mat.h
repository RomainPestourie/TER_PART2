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
  Eigen::VectorXd _Tvect;
  Eigen::MatrixXd _L1 , _L2 , _L3 , _L4 , _R;
  Eigen::MatrixXd _Rho , _A, _AStar , _Lambda , _RhoStar , _Xi , _XiStar;

  Eigen::VectorXd _FS , _FN , _FE , _FO, _Sm;

public:
Matrices(const double x_min, const double x_max, const double y_min, const double y_max,const int Nx, const int Ny);

void Flux();
void Rho(double t);
void RhoStar(double t);
void Xi();
void XiStar();
void Lambda();
void A();
void AStar();
void M();
void L1234();
void R();
void iteration();
void Sm();
void SaveSolPara(std::string nf);
Eigen::VectorXd GetT() {return _T;};
Eigen::SparseMatrix<double> GetM() {return _M;};
Eigen::MatrixXd GetRho() {return _Rho;};
Eigen::MatrixXd GetA() {return _A;};
Eigen::MatrixXd GetR() {return _R;};
Eigen::MatrixXd GetXi() {return _Xi;};
Eigen::MatrixXd GetLambda() {return _Lambda;};
Eigen::MatrixXd GetL1() {return _L1;};
Eigen::MatrixXd GetL2() {return _L2;};
Eigen::MatrixXd GetL3() {return _L3;};
Eigen::MatrixXd GetL4() {return _L4;};
Eigen::MatrixXd GetFlux() {return _f;};
};
