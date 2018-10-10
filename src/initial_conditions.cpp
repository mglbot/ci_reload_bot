#include "initial_conditions.hpp"
#include "quadrature.hpp"
#include <cmath>

void initial_conditions_1D(int const Lev_x, int const Lev_v, int const k,
                           PDE const pde, fk::matrix const FMWT_COMP_x,
                           int const FMWTxLDA, fk::matrix const FMWT_COMP_v,
                           int const FMWTvLDA, fk::vector &f_x2,
                           fk::vector &f_v2)
{
  // Get the Legendre-Gauss nodes (quad_x) and weights (quad_w) on the domain
  // [-1,+1] for performing quadrature.

  int const quad_num = 10;

  auto const [quad_x, quad_w] = []() -> std::array<fk::vector, 2> {
    fk::vector quad_x, quad_w;
    lgwt(quad_x, quad_w, quad_num, -1, 1);
    return {quad_x, quad_w};
  }();

  int const N = quad_x.size();

  // Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
  // to order k.

  fk::matrix const p_val = legendre(quad_x, k);

  // Get grid spacing for both  x and v.

  int const nx       = std::pow(2, Lev_x);
  double const hx    = pde.params.Lmax / nx;
  int const dof_1D_x = k * nx;
  fk::vector f_x(dof_1D_x);

  int const nv       = std::pow(2, Lev_v);
  double const hv    = 2.0 * pde.params.Vmax / nv;
  int const dof_1D_v = k * nv;
  fk::vector f_v(dof_1D_v);

  // Initial Condition for f_x and f_v

  fk::vector xi_x(N);
  fk::vector xi_v(N);
  for (auto i = 0; i <= nx - 1; ++i)
  {
    // Map quad_x from [-1,+1] to [0,LMax] physical domain.
    //
    for (auto j = 0; j < N; ++j)
    {
      xi_x(j) = hx * (quad_x(j) / 2.0 + 1.0 / 2.0 + i);
      xi_v(j) = hv * ((quad_x(j) + 1.0) / 2.0 + i) - pde.params.Vmax;
    }

    // Get the f(x) initial condition at the quadrature points.
    //
    fk::vector fxHere(N);
    fk::vector fvHere(N);
    for (auto j = 0; j < N; ++j)
    {
      fxHere(j) = pde.Fx_0(xi_x(j), pde.params); // FIXME shouldn't need params
      fvHere(j) = pde.Fv_0(xi_v(j), pde.params); // FIXME shouldn't need params
    }

    // Generate the coefficients for DG basis
    //
    for (auto thisk = 1; thisk <= k; ++thisk)
    {
      fk::vector this_k_legendre(N);
      fk::vector this_quad_x(N);
      fk::vector this_quad_v(N);

      for (auto m = 0; m < N; ++m)
      {
        int row            = thisk - 1;
        int col            = m;
        this_k_legendre(m) = p_val(row, col);

        this_quad_x(m) = (quad_w(m) * fxHere(m));
        this_quad_v(m) = (quad_w(m) * fvHere(m));
      }

      double tmp_x = this_k_legendre * this_quad_x;
      double tmp_v = this_k_legendre * this_quad_v;

      f_x(k * i + thisk - 1) = tmp_x * hx * std::sqrt(1.0 / hx) / 2.0;
      f_v(k * i + thisk - 1) = tmp_v * hv * std::sqrt(1.0 / hv) / 2.0;
    }
  }

  // Transfer to multi-DG bases

  // Matrix sizes
  // f_x (dof_1D_x,1)
  // f_v (dof_1D_v,1)
  // FMWTx (dof_1D_x,dof_1D_x)
  // FMWTv (dof_1D_v,dof_1D_v)

  f_x2.resize(dof_1D_x);
  f_v2.resize(dof_1D_v);

  f_v2 = f_v * FMWT_COMP_v;
  f_x2 = f_x * FMWT_COMP_x;
}
