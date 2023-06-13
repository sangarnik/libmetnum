#include <metnum.hpp>

void p1(metnum::ODESolverMethods method, bool plot = false) {
  metnum::ODESolverInput in;
  in.syst = {
    [](const metnum::MultiVarFuncArg& y) {return 3*y[0] - 2*y[1];},
    [](const metnum::MultiVarFuncArg& y) {return 5*y[0] - 4*y[1];}
  };
  in.init = {3, 6};
  in.h = 0.2;
  in.x_max = 2;

  std::cout << "Ejercicio 1): "; 

  if (method == metnum::EULER_SIMPLE) {
    std::cout << "Euler Simple" << std::endl;
    metnum::ODESolverEulerSimple out(in);
    out.calculate();
    out.print();
    if (plot)
      out.plot();
  }
  else if (method == metnum::RUNGE_KUTTA_4) {
    std::cout << "Runge Kutta 4" << std::endl;
    metnum::ODESolverRungeKutta4 out(in);
    out.calculate();
    out.print();
    if (plot)
      out.plot();
  }
}

void p2(metnum::ODESolverMethods method, bool plot = false) {
  metnum::ODESolverInput in;
  in.syst = {
    [](const metnum::MultiVarFuncArg& y) {return y[2];},
    [](const metnum::MultiVarFuncArg& y) {return 2*exp(y[0]) - 2*y[1] - y[2];}
  };
  in.h = 0.1;
  in.x_max = 2;
  in.init = {0, 1};

  std::cout << "Ejercicio 2): "; 

  if (method == metnum::EULER_SIMPLE) {
    std::cout << "Euler Simple" << std::endl;
    metnum::ODESolverEulerSimple out(in);
    out.calculate();
    out.print();
    if (plot)
      out.plot();
  }
  else if (method == metnum::RUNGE_KUTTA_4) {
    std::cout << "Runge Kutta 4" << std::endl;
    metnum::ODESolverRungeKutta4 out(in);
    out.calculate();
    out.print();
    if (plot)
      out.plot();
  }
}

int main() {
  //p1(metnum::EULER_SIMPLE, true);
  //p1(metnum::RUNGE_KUTTA_4, true);
  //p2(metnum::EULER_SIMPLE, true);
  p2(metnum::RUNGE_KUTTA_4, true);
}
