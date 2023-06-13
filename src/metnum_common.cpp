#include <metnum_common.hpp>

namespace metnum {
double deriv(const Function& f, double a, double delta) {
  return ( f(a+delta) - f(a) ) / delta;
}

double err(const double& a, const double& b) {
  return std::abs((b - a) / b );
}
}
