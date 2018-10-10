#ifndef PDE_HPP_
#define PDE_HPP_

#include <cmath>
#include <functional>

class PARAMS
{
public:
  double Lmin = 0.0;
  double Lmax = 20.0 * M_PI / 3.0;
  double Vmin = -13.0;
  double Vmax = +13.0;
};

class PDE
{
public:
  PARAMS params;

  unsigned int dim          = 2;
  unsigned int lev          = 3;
  unsigned int deg          = 2;
  unsigned int n_time_steps = 10;

  bool solvePoisson    = true;
  bool applySpecifiedE = false;
  bool checkAnalytic   = false;

  // Initial conditions
  //
  std::function<double(double, PARAMS)> Fx_0 = [&](double x, PARAMS p) {
    return 1 + 0.04 * std::cos(0.3 * x);
  };

  std::function<double(double, PARAMS)> Fv_0 = [&](double v, PARAMS p) {
    double np = 9.0 / (10.0 * std::sqrt(2.0 * M_PI));
    double nb = 2.0 / (10.0 * std::sqrt(2.0 * M_PI));
    double u  = 4.5;
    double vt = 0.5;

    return np * std::exp(-std::pow(v, 2) / 2.0) +
           nb * std::exp(-std::pow(v - u, 2) / (2.0 * std::pow(vt, 2)));
  };

  // Applied E field
  //
  std::function<double(double, PARAMS)> exactEx = [&](double x, PARAMS p) {
    return 0;
  };
  std::function<double(double, PARAMS)> exactEt = [&](double t, PARAMS p) {
    return 0;
  };

  // Exact solution
  //
  std::function<double(double, PARAMS)> exactFx = [&](double x, PARAMS p) {
    return x * 0;
  };
  std::function<double(double, PARAMS)> exactFv = [&](double v, PARAMS p) {
    return v * 0;
  };
  std::function<double(double, PARAMS)> exactFt = [&](double t, PARAMS p) {
    return t * 0;
  };

  // Source functions
  //
  std::function<double(double, PARAMS)> source1x = [&](double x, PARAMS p) {
    return x * 0;
  };
  std::function<double(double, PARAMS)> source1v = [&](double v, PARAMS p) {
    return v * 0;
  };
  std::function<double(double, PARAMS)> source1t = [&](double t, PARAMS p) {
    return t * 0;
  };

  std::function<double(double, PARAMS)> source2x = [&](double x, PARAMS p) {
    return x * 0;
  };
  std::function<double(double, PARAMS)> source2v = [&](double v, PARAMS p) {
    return v * 0;
  };
  std::function<double(double, PARAMS)> source2t = [&](double t, PARAMS p) {
    return t * 0;
  };

  std::function<double(double, PARAMS)> source3x = [&](double x, PARAMS p) {
    return x * 0;
  };
  std::function<double(double, PARAMS)> source3v = [&](double v, PARAMS p) {
    return v * 0;
  };
  std::function<double(double, PARAMS)> source3t = [&](double t, PARAMS p) {
    return t * 0;
  };
};
#endif
