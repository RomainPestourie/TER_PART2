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


class resol{
private:


public:

  void Newton();
  void Resolution();
  void Calc_Rho_etoile();
  void Calc_Rho_it();
  void DirectSolver();

};
