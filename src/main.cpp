
#include "operator_two_scale.hpp"
#include "pde.hpp"

#include <iostream>

int main(int argc, char **argv)
{
  PDE pde;
  auto params = pde.params;

  // Problem parameters

  std::cout << "--run begin--" << std::endl;

  if (argc == 3)
  {
    pde.deg = atoi(argv[1]);
    pde.lev = atoi(argv[2]);
    std::cout << "using supplied lev and degree" << std::endl;
  }
  else if (argc == 4)
  {
    pde.deg          = atoi(argv[1]);
    pde.lev          = atoi(argv[2]);
    pde.n_time_steps = atoi(argv[3]);
    std::cout << "using supplied lev, degree and n_time_steps" << std::endl;
  }
  else
  {
    std::cout << "using default lev and degree" << std::endl;
  }

  // int dim = pde.dim;
  int lev = pde.lev;
  int deg = pde.deg;
  // unsigned long n_time_steps = pde.n_time_steps;

  std::cout << "deg :" << deg << ": " << std::endl;
  std::cout << "lev :" << lev << ": " << std::endl;

  int LevX = lev;
  int LevV = lev;

  // Time step
  //
  double dt =
      params.Lmax / std::pow(2.0, LevX) / params.Vmax / (2.0 * deg + 1.0);

  std::cout << "dt :" << dt << ": " << std::endl;

  // Step 1.1. Set up matrices for multi-wavelet
  //
  fk::matrix FMWT_COMP_x_2d =
      operator_two_scale(deg, static_cast<int>(pow(2, LevX)));
  fk::matrix FMWT_COMP_v_2d =
      operator_two_scale(deg, static_cast<int>(pow(2, LevV)));
}
