#include <metnum_diffeq.hpp>

namespace metnum {
ODESolverBase::ODESolverBase(const ODESolverInput& in) {
  /* Init */
  // Load initial values
  this->values.push_back({0});
  for (int i = 0; i < in.init.size(); ++i)
    this->values.push_back({in.init[i]});

  // Load the rest of inputs
  this->syst = in.syst;
  this->h = in.h;
  this->x_max = in.x_max;
  
  std::cout << "Inited" << std::endl;
  //calculate(); -> constructor cannot be overriden, dont do this
}

void ODESolverBase::print() {
  /* Print calculated values */
  std::cout << "Printing Values:" << std::endl;
  for (int i = 0; i <= getRange(); ++i) {
    for (int j = 1; j < getEqNum()+1; ++j) {
      std::cout << "y_" << j << "(" << values[0][i] << ") = " << values[j][i];
      if (j != getEqNum())
        std::cout << " ; ";
    }
    std::cout << std::endl;
  }
}

void ODESolverBase::calculate() {
  /* Calculate values */
  // y_{i+1} = y_i + \phi.h -> \phi depende del metodo
  // x_{i+1} = x_i + h -> \phi = 1 siempre
  for (int i = 0; i < getRange(); ++i) {
    int eq = 0;
    for (auto& phi : getIncrement()) {
      values[eq].push_back(values[eq].back() + (phi*h) );
      ++eq;
    }
  }
  std::cout << "Calculated" << std::endl;
}

std::vector<double> ODESolverBase::getLastVals() {
  /* Return latest values before iteration */
  std::vector<double> arg; // arg[0] = x, arg[n] = y_n

  for (int i = 0; i <= getEqNum(); ++i)
    arg.push_back(values[i].back());

  return arg;
};

void ODESolverBase::writeToFile(std::string filename) {
  /* Write values to file */
  std::filesystem::path f {filename};
  if (std::filesystem::exists(f)) // Update old file if exists
    std::remove(filename.c_str());

  std::ofstream outFile(filename);

  for (int i = 0; i <= getRange(); ++i) {
    for (int j = 0; j <= getEqNum(); ++j) {
      outFile << values[j][i];
      if (j != getEqNum())
        outFile << "\t";
      else
        outFile << "\n";
    }
  }

  outFile.close();
}

void ODESolverBase::plot() {
  /* Graph values with gnuplot */
  GnuplotPipe gplot;
  std::string plotString = "plot ";
  writeToFile("plot.txt"); // Load values into file

  gplot.sendLine("set xzeroaxis");
  gplot.sendLine("set yzeroaxis");
  gplot.sendLine("set grid");
  gplot.sendLine("set terminal qt size 959, 650");
  
  for (int i = 0; i < getEqNum(); ++i) {
    plotString += "'plot.txt' u 0:"+std::to_string(i+2)+":xtic(1) w l title 'y_"+std::to_string(i+1)+"'";
    if (i != getEqNum()-1)
      plotString += ",";
  }
  
  gplot.sendLine(plotString);
}

void ODESolverBase::comparePlots(std::vector<ODESolverBase> plots, std::string filename) { 
  /* Write N solutions to file and graph with gnuplot */
  int eqNumber, eqRange;
  ODEValues out;
  GnuplotPipe gplot;
  std::vector<ODEValues> vals;
  std::string plotString = "plot ";

  std::filesystem::path f {filename};
  if (std::filesystem::exists(f))
    std::remove(filename.c_str());

  std::ofstream outFile(filename);

  gplot.sendLine("set xzeroaxis");
  gplot.sendLine("set yzeroaxis");
  gplot.sendLine("set grid");
  gplot.sendLine("set terminal qt size 959, 650");

  for (auto& plot : plots)
    vals.push_back(plot.getValues());
  out = vals.front();
  
  for (int i = 1; i < vals.size(); ++i) {
    vals[i].erase(vals[i].begin());
    for (int j = 0; j < vals[i].size(); ++j)
      out.push_back(vals[i][j]);
  }

  eqNumber = out.size();
  eqRange = plots.front().getRange();
  
  for (int i = 0; i <= eqRange; ++i)
    for (int j = 0; j < eqNumber; ++j) {
      outFile <<  out[j][i];
      if (j != eqNumber-1)
        outFile << "\t";
      else
        outFile << "\n";
    }
  outFile.close();
  
  for (int i = 0; i < eqNumber-1; ++i) {
    plotString += "'"+filename+"' u 0:"+std::to_string(i+2)+":xtic(1) w l title 'y_"+std::to_string(i+1)+"'";
    if (i != eqNumber-1)
      plotString += ",";
  }
  gplot.sendLine(plotString);
}

IncrVector ODESolverEulerSimple::getIncrement() {
  /* Calculate \phi for Euler Simpl */
  MultiVarFuncArg arg = getLastVals();
  IncrVector out(arg.size());

  // \phi(x,y,h) = f(x_i, y_i)
  // \phi(x,y_1,...,y_n,h) = f(x_i, y_1_i,...,y_n_i)
  out[0] = 1; // x_i incrementa en 1*h
  for (int i = 0; i < getEqNum(); ++i)
    // y_i incrementan en \phi*h
    out[i+1] = syst[i](arg);
  
  return out;
}

IncrVector ODESolverRungeKutta4::getIncrement() {
  /* Calculate \phi for Runge Kutta */
  std::vector<double> y = getLastVals(); // y[0] = x, y[n] = y_n
  MultiVarFuncArg arg = y;
  IncrVector out(arg.size());
  KValues k(getEqNum());

  // k_1 = f(x_i, y_i)
  for (int i = 0; i < getEqNum(); ++i)
    k[i].push_back(syst[i](arg));

  // k_2 = f( x_i + h/2, y_i + (k_1*h*0.5) )
  arg[0] = y[0] + (h * 0.5);
  for (int i = 1; i < arg.size(); ++i)
    arg[i] = y[i] + (k[i-1][0] * h * 0.5);
  for (int i = 0; i < getEqNum(); ++i)
    k[i].push_back(syst[i](arg));

  // k_3 = f( x_i + h/2, y_i + (k_2*h*0.5) )
  //arg[0] = y[0] + (h * 0.5); <- arg[0] no cambia
  for (int i = 1; i < arg.size(); ++i)
    arg[i] = y[i] + (k[i-1][1] * h * 0.5);
  for (int i = 0; i < getEqNum(); ++i)
    k[i].push_back(syst[i](arg));
  
  // k_4 = f( x_i + h, y_i + (k_2*h) )
  arg[0] = y[0] + h;
  for (int i = 1; i < arg.size(); ++i)
    arg[i] = y[i] + (k[i-1][2] * h);
  for (int i = 0; i < getEqNum(); ++i)
    k[i].push_back(syst[i](arg));

  // \phi_i = (1/6)*(k_i,1 + 2*k_i,2 + 2*k_i,3 + k_i,4)
  out[0] = 1; // x_i incrementa en 1*h
  for (int i = 0; i < getEqNum(); ++i)
    // y_i incrementan en \phi*h
    out[i+1] = (k[i][0] + 2*(k[i][1] + k[i][2]) + k[i][3])/6.0;

  return out;
}
}
