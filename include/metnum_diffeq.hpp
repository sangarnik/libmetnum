#pragma once

#include <fstream>
#include <filesystem>

#include <metnum_common.hpp>
#include <metnum_gnuplot.hpp>

namespace metnum {
typedef std::vector<double> ODEInitValues, IncrVector;
typedef std::vector<std::vector<double>> ODEValues, KValues;
typedef std::vector<MultiVarFunc> ODESystem;

enum ODESolverMethods {
  EULER_SIMPLE,
  RUNGE_KUTTA_4
};

/* Input for ODE Solvers */
struct ODESolverInput {
  ODEInitValues init;
  ODESystem syst;
  double h;
  double x_max;
};

/* Base class for ODE Solvers */
class ODESolverBase {
public:
  /* Init */
  ODESolverBase(const ODESolverInput& in);
  /* Write N solutions to file and graph with gnuplot */
  static void comparePlots(std::vector<ODESolverBase> plots, std::string filename = "plotcompare.txt");
  /* Get total number of equations */
  int getEqNum() {return this->syst.size();};
  /* Get calculated values */
  ODEValues getValues() {return this->values;};
  /* Get solution range in the x axis */
  double getRange() {return (this->x_max/this->h);};
  /* Print calculated values */
  void print();
  /* Graph values with gnuplot */
  void plot();
  /* Calculate values */
  void calculate();
  /* Write values to file */
  void writeToFile(std::string filename = "plot.txt");
protected:
  double h;
  double x_max;
  ODESystem syst;
  ODEValues values;
  /* Overridable function to calculate \phi (increment) */
  virtual IncrVector getIncrement() {return IncrVector();};
  /* Return latest values before iteration */
  std::vector<double> getLastVals();
};

/* Euler Simple ODE solver method */
class ODESolverEulerSimple : public ODESolverBase {
public:
  ODESolverEulerSimple(const ODESolverInput& in) : ODESolverBase(in) {};
protected:
  /* Calculate \phi for Euler Simple */
  IncrVector getIncrement() override;
};

/* Runge Kutta 4 ODE solver method */
class ODESolverRungeKutta4 : public ODESolverBase {
public:
  ODESolverRungeKutta4(const ODESolverInput& in) : ODESolverBase(in) {};
protected:
  /* Calculate \phi for Runge Kutta */
  IncrVector getIncrement() override;
};
}
