#pragma once

#include <metnum_common.hpp>

namespace metnum {
/* Input for closed root finder methods */
struct ClosedRootInput {
  double a;
  double b;
  double err;
  int iter;
};
/* Input for open root finder methods */
struct OpenRootInput {
  double x0;
  double err;
  int iter;
};

enum ClosedRootFinderMethods {
  BISECCION,
  REGULA_FALSI,
  PUNTO_FIJO
};

enum OpenRootFinderMethods {
  NEWTON,
  SECANTE
};

/* Calculate optimal iterations for root finders */
int rootStopIter(double a, double b, double err);
/* Raíz por bisección (Closed) */
int rootBisecc(const Function& f, const ClosedRootInput& in, double& out);
/* Raíz por regula falsi (Closed) */
int rootRegFal(const Function& f, const ClosedRootInput& in, double& out);
/* Raiz por punto fijo (despejar ecuación con respecto a x) (Open) */
int rootPtFijo(const Equation& f, const OpenRootInput& in, double& out);
/* Raíz por método de Newton (Open) */
int rootNewton(const Function& f, const OpenRootInput& in, double& out);
/* Raíz por método de la secante (Open) */
int rootSecant(const Function& f, const OpenRootInput& in, double& out);
}
