#pragma once

#include <armadillo>

#include <metnum_common.hpp>

namespace metnum {
typedef std::vector<double> VecDouble;
typedef std::vector<VecDouble> MatDouble;

struct EquationMat {
  MatDouble coef;
  VecDouble res;
};

struct IterMetInput {
  EquationMat mat;
  int iter;
  double err;
  VecDouble x_0;
};

enum DirectEqSolverMethods {
  GAUSS_SIMPLE,
  GAUSS_JORDAN
};

enum IterEqSolverMethods {
  JACOBI,
  GAUSS_SEIDEL
};

enum MatCondic {
  COND_OK,
  COND_DET_ZERO,
  COND_BIG_INV,
  COND_BIG_IDENT
};

void systGaussSim(const EquationMat& mat, VecDouble& out);
void systGaussJor(const EquationMat& mat, VecDouble& out);
int systJacobi(const IterMetInput& in, VecDouble& out);
int systGaussSei(const IterMetInput& in, VecDouble& out);
}

namespace matrixutil {
arma::mat toArmaMat(const metnum::MatDouble& mat);
void printVec(const metnum::VecDouble& vec);
void printEqMat(const metnum::EquationMat& mat);
double findMax(const metnum::MatDouble& mat);
double findMax(const metnum::VecDouble& vec);
double det(const metnum::MatDouble& mat);
metnum::EquationMat normalise(const metnum::EquationMat& mat);
metnum::MatCondic checkCond(const metnum::EquationMat& mat);
void printCond(const metnum::EquationMat& mat);
bool willConverge(const metnum::EquationMat& mat);
//double getErr(const double& a, const double& b);
}
