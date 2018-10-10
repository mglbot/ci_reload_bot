
#include "quadrature.hpp"

#define _USE_MATH_DEFINES
#include "matlab_utilities.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>

#include <iostream>
//-----------------------------------------------------------------------------
//
// initial c++ implementation by Tyler McDaniel
//
// output vectors (x and w) are allocated by caller and passed by ref
//
//
// --- (matlab description) ---
//
// function [x,w]=lgwt(N,a,b)
//  lgwt.m
//
//  This script is for computing definite integrals using Legendre-Gauss
//  Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
//  [a,b] with truncation order N
//
//  Suppose you have a continuous function f(x) which is defined on [a,b]
//  which you can evaluate at any x in [a,b]. Simply evaluate it at all of
//  the values contained in the x vector to obtain a vector f. Then compute
//  the definite integral using sum(f.*w);
//
//  Written by Greg von Winckel - 02/25/2004
//
//-----------------------------------------------------------------------------

void lgwt(fk::vector &x_out, fk::vector &w_out, unsigned int const n,
          int const a, int const b)
{
  assert(n > 0);

  // prepare out vectors
  std::vector<double> x, w;
  x.resize(n);
  x_out.resize(n);
  w.resize(n);
  w_out.resize(n);

  // N=N-1;
  // N1=N+1; N2=N+2;
  unsigned int const N  = n - 1;
  unsigned int const N1 = N + 1;
  unsigned int const N2 = N + 2;

  // xu=linspace(-1,1,N1)';
  std::vector<double> xu = linspace(-1.0, 1.0, N1);

  //% Initial guess
  // y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
  std::vector<double> y = linspace(0, N, N1);
  std::transform(y.begin(), y.end(), y.begin(), [&](double &elem) {
    return std::cos((2 * elem + 1) * M_PI / (2 * N + 2));
  });

  std::vector<double> y2(xu);
  std::transform(y2.begin(), y2.end(), y2.begin(), [&](double &elem) {
    return (0.27 / N1) * std::sin(M_PI * elem * N / N2);
  });

  std::transform(y.begin(), y.end(), y2.begin(), y.begin(),
                 std::plus<double>());

  //% Legendre-Gauss Vandermonde Matrix
  // L=zeros(N1,N2);
  std::vector<double> L(N1 * N2, 0.0);

  //% Derivative of LGVM
  // Lp=zeros(N1,N2);
  // resized, size in MATLAB is in error -TM
  std::vector<double> Lp(N1, 0.0);

  //% Compute the zeros of the N+1 Legendre Polynomial
  // y0=2
  std::vector<double> y0(N1, 2.0);
  double eps = std::numeric_limits<double>::epsilon();

  //% Iterate until new points are uniformly within epsilon of old points
  // while max(abs(y-y0))>eps
  std::vector<double> diff(N1);
  std::transform(y.begin(), y.end(), y0.begin(), diff.begin(),
                 [&](double &y_elem, double &y0_elem) {
                   return std::fabs(y_elem - y0_elem);
                 });

  while (*std::max_element(diff.begin(), diff.end()) > eps)
  {
    // L(:,1)=1;
    // Lp(:,1)=0;
    // L(:,2)=y;
    // Lp(:,2)=1;
    for (unsigned int i = 0; i < N1; ++i)
    {
      L[i]      = 1.0;
      L[i + N1] = y[i];
      // Lp column one initialization performed above, second column
      // initalization is unnecessary
    }

    // for k=2:N1
    // iterating over the remaining columns of L
    for (unsigned int i = 1; i < N1; ++i)
    {
      // calc column offsets -TM
      auto current  = i * N1;
      auto previous = current - N1;
      auto next     = current + N1;
      for (unsigned int j = 0; j < N1; ++j)
      {
        // L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
        L[next + j] =
            ((2 * (i + 1) - 1) * y[j] * L[current + j] - i * L[previous + j]) /
            (i + 1);
      }
    }
    // end

    // Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);
    for (unsigned int i = 0; i < N1; ++i)
    {
      Lp[i] = N2 * (L[(N1 - 1) * N1 + i] - y[i] * L[(N2 - 1) * N1 + i]) /
              (1 - pow(y[i], 2));
    }

    // y0=y;
    y0 = y;

    // y=y0-L(:,N2)./Lp;
    for (unsigned int i = 0; i < N1; ++i)
    {
      y[i] = y0[i] - L[(N2 - 1) * N1 + i] / Lp[i];
    }

    std::transform(y.begin(), y.end(), y0.begin(), diff.begin(),
                   [](double &y_elem, double &y0_elem) {
                     return std::fabs(y_elem - y0_elem);
                   });
    // end
  }

  //% Linear map from[-1,1] to [a,b]
  // x=(a*(1-y)+b*(1+y))/2;
  std::transform(y.begin(), y.end(), x.begin(), [&](double &elem) {
    return (a * (1 - elem) + b * (1 + elem)) / 2;
  });

  //% Compute the weights
  // w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
  std::transform(y.begin(), y.end(), Lp.begin(), w.begin(),
                 [&](double &y_elem, double &Lp_elem) {
                   return (b - a) / ((1 - pow(y_elem, 2)) * pow(Lp_elem, 2)) *
                          pow(static_cast<double>(N2) / N1, 2);
                 });

  // x=x(end:-1:1);
  // w=w(end:-1:1);
  std::reverse(x.begin(), x.end());
  std::reverse(w.begin(), w.end());

  x_out = x;
  w_out = w;
}

//% Legendre Polynomials with degree k on [-1,1]
fk::matrix legendre(fk::vector const x_in, unsigned int const k)
{
  // legendre: currently harcoded to support degree 2 - 7"
  assert(k <= 7);
  assert(k >= 2);

  std::vector<double> x = x_in.to_std();
  assert(!x.empty());

  std::vector<double> v(x.size());
  fk::matrix m(k, x.size());

  switch (k)
  {
  case 7:
    // v(:,7)=1/16*(231*x.^6-315*x.^4+105*x.^2-5)*sqrt(7*2-1);
    std::transform(x.begin(), x.end(), v.begin(), [](double &elem) {
      return 0.0625 *
             (231 * pow(elem, 6) - 315 * pow(elem, 4) + 105 * pow(elem, 2) -
              5) *
             sqrt(13);
    });
    m.update_row(6, v);
  case 6:
    // v(:,6)=1/8*(63*x.^5-70*x.^3+15*x)*sqrt(6*2-1);
    std::transform(x.begin(), x.end(), v.begin(), [](double &elem) {
      return 0.125 * (63 * pow(elem, 5) - 70 * pow(elem, 3) + 15 * elem) *
             sqrt(11);
    });
    m.update_row(5, v);
  case 5:
    // v(:,5)=1/8*(35*x.^4-30*x.^2+3)*sqrt(5*2-1);
    std::transform(x.begin(), x.end(), v.begin(), [](double &elem) {
      return 0.125 * (35 * pow(elem, 4) - 30 * pow(elem, 2) + 3) * sqrt(9);
    });
    m.update_row(4, v);
  case 4:
    // v(:,4)=1/2*(5*x.^3-3*x)*sqrt(4*2-1);
    std::transform(x.begin(), x.end(), v.begin(), [](double &elem) {
      return 0.5 * (5 * pow(elem, 3) - 3 * elem) * sqrt(7);
    });
    m.update_row(3, v);
  case 3:
    // v(:,3)=1/2*(3*x.^2-1)*sqrt(3*2-1);
    std::transform(x.begin(), x.end(), v.begin(), [](double &elem) {
      return 0.5 * (3 * pow(elem, 2) - 1) * sqrt(5);
    });
    m.update_row(2, v);
  case 2:
    // v(:,2)=x*sqrt(2^2-1);
    std::transform(x.begin(), x.end(), v.begin(),
                   [](double &elem) { return elem * sqrt(3); });
    m.update_row(1, v);
    // v(:,1)=x-x+1;
    v = std::vector<double>(x.size(), 1.0);
    m.update_row(0, v);
  }

  return m;
}

//% Legendre Polynomials with degree k on [-1,1]
fk::matrix dlegendre(fk::vector const x_in, unsigned int const k)
{
  // legendre: currently harcoded to support degree 2 - 7"
  assert(k <= 7);
  assert(k >= 2);

  std::vector<double> x = x_in.to_std();
  assert(!x.empty());

  std::vector<double> v(x.size());
  fk::matrix m(k, x.size());

  switch (k)
  {
  case 7:
    // v(:,7)=1/16*(231*6*x.^5-315*4*x.^3+105*2*x)*sqrt(7*2-1);
    std::transform(x.begin(), x.end(), v.begin(), [](double &elem) {
      return 0.0625 *
             (231 * 6 * pow(elem, 5) - 315 * 4 * pow(elem, 3) + 105 + 2 +
              elem) *
             sqrt(13);
    });
    m.update_row(6, v);
  case 6:
    // v(:,6)=1/8*(63*5*x.^4-70*3*x.^2+15)*sqrt(6*2-1);
    std::transform(x.begin(), x.end(), v.begin(), [](double &elem) {
      return 0.125 * (63 * 5 * pow(elem, 4) - 70 * 3 * pow(elem, 2) + 15) *
             sqrt(11);
    });
    m.update_row(5, v);
  case 5:
    // v(:,5)=1/8*(35*4*x.^3-30*2*x)*sqrt(5*2-1);
    std::transform(x.begin(), x.end(), v.begin(), [](double &elem) {
      return 0.125 * (35 * 4 * pow(elem, 3) - 30 * 2 * elem) * sqrt(9);
    });
    m.update_row(4, v);
  case 4:
    // v(:,4)=1/2*(5*3*x.^2-3)*sqrt(4*2-1);
    std::transform(x.begin(), x.end(), v.begin(), [](double &elem) {
      return 0.5 * (15 * pow(elem, 2) - 3) * sqrt(7);
    });
    m.update_row(3, v);
  case 3:
    // v(:,3)=1/2*(3*2*x)*sqrt(3*2-1);
    std::transform(x.begin(), x.end(), v.begin(),
                   [](double &elem) { return 0.5 * (6 * elem) * sqrt(5); });
    m.update_row(2, v);
  case 2:
    // v(:,2)=x-x+1*sqrt(2*2-1);
    std::transform(x.begin(), x.end(), v.begin(),
                   [](double &elem) { return sqrt(3); });
    m.update_row(1, v);
    // v(:,1)=x-x;
    v = std::vector<double>(x.size(), 0.0);
    m.update_row(0, v);
  }

  return m;
}
