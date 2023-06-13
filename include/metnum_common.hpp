#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <functional>


namespace metnum {
typedef std::function<double(double)> Function, Equation;
typedef std::vector<double> MultiVarFuncArg;
typedef std::function<double(MultiVarFuncArg)> MultiVarFunc;

/* Get derivative for the function f(x) in x=a */
double deriv(const Function& f, double a, double delta = 0.001);
/* Get relative error between a and b */
double err(const double& a, const double& b);
}
