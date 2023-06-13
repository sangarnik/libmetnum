#include <metnum_rootfinder.hpp>

namespace metnum {
int rootStopIter(double a, double b, double err) {
  return round( log( (b-a)/(err) / log(2) ) );
}

//int getRootBiseccion(double* out, double (*f)(double), ClosedRootInput& in) {
int rootBisecc(const Function& f, const ClosedRootInput& in, double& out) {
  int iter = in.iter;
  double currerr = 1;
  double a = in.a;
  double b = in.b;
  double c;

  if ( f(a)*f(b) < 0 ) {
    for ( int i = 0; i < iter; ++i ) {
      double cprev = c;
      c = (a+b)*0.5f;
      //std::cout << c << std::endl;
      if ( f(a)*f(c) < 0 )
        b = c;
      else
        a = c;

      if ( f(c) == 0 || currerr < in.err )
        break;
      else
        currerr = metnum::err(cprev, c);
    }
    out = c;
    return 0;
  }
  else
    return 1;
}

//int getRootRegulaFalsi(double* out, double (*f)(double), ClosedRootInput& in) {
int rootRegFal(const Function& f, const ClosedRootInput& in, double& out) {
  int iter = in.iter;
  double currerr = 1;
  double a = in.a;
  double b = in.b;
  double c;

  if ( f(a)*f(b) < 0 ) {
    for ( int i = 0; i < iter; ++i ) {
      double cprev = c;
      c = b - ( ( f(b)*(a-b) )/( f(a)-f(b) ) );
      //std::cout << c << std::endl;
      if ( f(a)*f(c) < 0 )
        b = c;
      else
        a = c;

      if ( f(c) == 0 || currerr < in.err )
        break;
      else
        currerr = metnum::err(cprev, c);
    }
    out = c;
    return 0;
  }
  else
    return 1;
}

//int getRootPuntoFijo(double* out, double (*f)(double), OpenRootInput& in) {
int rootPtFijo(const Equation& f, const OpenRootInput& in, double& out) {
  int iter = in.iter;
  double currerr = 1;
  double xr = in.x0;
  //std::cout << deriv(f,xr) << std::endl;

  if (std::abs(deriv(f, xr)) > 1)
    return 1;

  for ( int i = 0; i < iter; ++i ) {
    double x_i = xr;
    xr = f(x_i);
    //std::cout << xr << std::endl;
    if ( xr == 0 || currerr < in.err)
      break;
    else
      currerr = metnum::err(x_i, xr);
  }
  out = xr;

  return 0;
}

//int getRootMetNewton(double* out, double(*f)(double), OpenRootInput& in) {
int rootNewton(const Function& f, const OpenRootInput& in, double& out) {
  int iter = in.iter;
  double currerr = 1;
  double xr = in.x0;
  
  for ( int i = 0; i < iter; ++i ) {
    double x_i = xr;
    xr = x_i - (f(xr) / deriv(f, xr));
    //std::cout << xr << std::endl;
    if ( xr == 0 || currerr < in.err )
      break;
    else
      currerr = metnum::err(x_i, xr);
  }
  out = xr;

  return 0;
}

//int getRootMetSecante(double* out, double (*f)(double), OpenRootInput& in) {
int rootSecant(const Function& f, const OpenRootInput& in, double& out) {
  int iter = in.iter;
  double currerr = 1;
  double xr = in.x0;
  double x_im1 = xr+0.1;
  
  for ( int i = 0; i < iter-1; ++i ) {
    double x_i = xr;
    xr = x_i - ((f(x_i)*(x_im1-x_i))/(f(x_im1)-f(x_i)));
    x_im1 = x_i;
    //std::cout << xr << std::endl;
    if ( xr == 0 || currerr < in.err )
      break;
    else
      currerr = metnum::err(x_i, xr);
  }
  out = xr;

  return 0;
}
}
