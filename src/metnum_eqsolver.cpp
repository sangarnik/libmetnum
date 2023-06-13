#include <metnum_eqsolver.hpp>

namespace metnum {
void systGaussSim(const EquationMat& mat, VecDouble& out) {
  MatDouble a = mat.coef;
  VecDouble b = mat.res;
  VecDouble x = VecDouble(mat.coef.size());
  std::fill(x.begin(), x.end(), 0);
  int n = b.size();

  // Escalonar matriz
  for (int k = 0; k < n-1; ++k)
    for (int i = k+1; i < n; ++i) {
      double factor = a[i][k] / a[k][k];
      for (int j = k+1; j < n; ++j)
        a[i][j] -= (factor*a[k][j]);
      a[i][k] = 0;
      b[i] -= (factor*b[k]);
    }

  // Sustitución hacia atras
  x[n-1] = b[n-1]/a[n-1][n-1];
  for (int i = n-1; i >= 0; --i) {
    double sum = b[i];
    for (int j = i+1; j < n; ++j)
      sum -= a[i][j]*x[j];
    x[i] = sum / a[i][i];
  }

  out = x;
}

void systGaussJor(const EquationMat& mat, VecDouble& out) {
  MatDouble a = mat.coef;
  VecDouble b = mat.res;
  int n = b.size();
  
  for (int k = 0; k < n; ++k) {
    double factor = a[k][k];
    // Normalizar fila
    std::for_each(a[k].begin()+k, a[k].end(), [factor](double& x){ x/= factor;});
    b[k] /= factor;

    for (int i = 0; i < n; ++i) {
      // Reducir filas restantes
      if (i == k)
        continue;
      
      double sfactor = a[i][k];
      b[i] -= b[k]*sfactor;
      for (int j = k; j < n; ++j)
        a[i][j] -= a[k][j]*sfactor;
    }
  }

  out = b; 
}

//int systEqJacobi(const EquationMat& mat, VecDouble& out, int iter, double err, VecDouble x_0) {
int systJacobi(const IterMetInput& in, VecDouble& out) {
  int n = in.x_0.size();
  int finaliter = in.iter;
  MatDouble a = in.mat.coef;
  VecDouble b = in.mat.res;
  VecDouble x_0 = in.x_0;
  VecDouble x_1 = x_0;

  for (int k = 0; k < in.iter; ++k) {
    for (int i = 0; i < n; ++i) {
      double sum = 0;
      for (int j = 0; j < n; ++j) {
        if (j == i)
          continue;
        sum += a[i][j]*x_0[j];
      }
      x_1[i] = (b[i] - sum) / a[i][i];
    }
    if (metnum::err(x_0[0], x_1[0]) < in.err) {
      //std::cout << k << std::endl;
      finaliter = k;
      break;
    }
    x_0 = x_1;
  }

  out = x_1;
  return finaliter;
}

//int systEqGaussSeidel(const EquationMat& mat, VecDouble& out, int iter, double err, VecDouble x_0) {
int systGaussSei(const IterMetInput& in, VecDouble& out) {
  int n = in.x_0.size();
  int finaliter = in.iter;
  MatDouble a = in.mat.coef;
  VecDouble b = in.mat.res;
  VecDouble x_0 = in.x_0;

  for (int k = 0; k < in.iter; ++k) {
    double old = x_0[0];
    for (int i = 0; i < n; ++i) {
      double sum = 0;
      for (int j = 0; j < n; ++j) {
        if (j == i)
          continue;
        sum += a[i][j]*x_0[j];
      }
      x_0[i] = (b[i] - sum) / a[i][i];
    }
    if (metnum::err(old, x_0[0]) < in.err) {
      //std::cout << k << std::endl;
      finaliter = k;
      break;
    }
  }

  out = x_0;
  return finaliter;
}
}

namespace matrixutil {
arma::mat toArmaMat(const metnum::MatDouble& mat) {
  int n = mat.size();
  int m = mat[0].size();
  arma::mat out(n,m);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      out(i,j) = mat[i][j];
    }
  }
  
  return out;
}

void printVec(const metnum::VecDouble& vec) {
  for (auto& item : vec)
    std::cout << item << std::endl;
}

void printEqMat(const metnum::EquationMat& mat) {
  std::cout << toArmaMat(mat.coef);
  printVec(mat.res);

  /*
  for (int i = 0; i < item.size(); ++i) {
    std::string out = "";
    for (int j = 0; j < coefs.size(); ++j) {
      out += std::to_string(coefs[i][j]) + "*x_" + std::to_string(j);
      if ( j != coefs.size() -1)
        out += " + ";
      else
        out += " = ";
    }
    std::cout << out << item[i] << std::endl;
  }
  std::cout << std::endl;
 */
}

double findMax(const metnum::MatDouble& mat ) {
  return toArmaMat(mat).max();
}

double findMax(const metnum::VecDouble& vec) {
  return arma::max(arma::vec(vec));
}

double det(const metnum::MatDouble& mat) {
  return arma::det(toArmaMat(mat));
}

metnum::EquationMat normalise(const metnum::EquationMat& mat) {
  double max = findMax(mat.coef);
  metnum::EquationMat out;
  out.coef = mat.coef;
  out.res = mat.res;

  for (auto& row : out.coef) {
    for (auto& item : row) {
      item = item/max;
    }
  }

  for (auto& item : out.res) {
    item = item/max;
  }
  return out;
}

/*
 * Corroborar mal condicionamiento:
 * 1) Normalizar ecuaciones (a partir del coeficiente mas grande)
 * 2) Det != 0
 * 3) Si existen elementos de la matriz inversa muy mayores a 1, posible mal condicionamiento [A]^-1 >>> 1
 * 4) Multiplicar A^-1 . A ~= I si no es cercano a la matriz identidad, mal condicionamiento
*/

metnum::MatCondic checkCond(const metnum::EquationMat& mat) {
  int mat_inv_treshold = 10;
  int mult_treshold = 2;
  arma::mat arma_mat = toArmaMat(mat.coef);
  arma::mat arma_mat_inv = arma::inv(arma_mat);
  arma::mat mult = arma_mat * arma_mat_inv;

  if (det(mat.coef) == 0)
    return metnum::COND_DET_ZERO;

  std::cout << arma_mat_inv;
  for (auto& item : arma_mat_inv) {
    if (item > mat_inv_treshold)
      return metnum::COND_BIG_INV;
  }

  for (auto& item : mult) {
    if (item > mult_treshold)
      return metnum::COND_BIG_IDENT;
  }

  return metnum::COND_OK;
}

void printCond(const metnum::EquationMat& mat) {
   switch (checkCond(mat)) {
    case metnum::COND_BIG_INV:
      std::cout << "COND_BIG_INV" << std::endl;
      break;
    case metnum::COND_DET_ZERO:
      std::cout << "COND_DET_ZERO" << std::endl;
      break;
    case metnum::COND_BIG_IDENT:
      std::cout << "COND_BIG_IDENT" << std::endl;
    default:
      std::cout << "COND_OK" << std::endl;
  }
}

/* 
* Condición convergencia seidel y jacobi
* |a_ii| > |\sum_{j=1 j != i}^{n} a_{ij}|
*/

bool willConverge(const metnum::EquationMat& mat) {
  bool flag = true;

  for (int i = 0; i < mat.coef.size(); ++i) {
    double curr = mat.coef[i][i];
    double sum = 0;
    for (int j = 0; j < mat.coef[0].size(); ++j) {
      if (j == i)
        continue;
      sum += std::abs(mat.coef[i][j]);
    }
    if (std::abs(curr) <= sum) {
      flag = false;
      break;
    }
  }

  return flag;
}

/*
double err(const double& a, const double& b){
  return (std::abs((b-a)/b));
}
*/
}
