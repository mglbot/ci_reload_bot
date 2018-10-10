#ifndef INITIAL_CONDITIONS_HPP_
#define INITIAL_CONDITIONS_HPP_

#include "pde.hpp"
#include "tensors.hpp"

void initial_conditions_1D(int const Lev_x, int const Lev_v, int const k,
                           PDE const pde, fk::matrix const FMWT_COMP_x,
                           int const FMWTxLDA, fk::matrix const FMWT_COMP_v,
                           int const FMWTvLDA, fk::vector &f_x2,
                           fk::vector &f_v2);

#endif
