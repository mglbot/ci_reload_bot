

//-----------------------------------------------------------------------------
//
// quadrature.hpp
//
//-----------------------------------------------------------------------------

#ifndef _quadrature_h_
#define _quadrature_h_

#include "tensors.hpp"
#include <vector>

void lgwt(fk::vector &x, fk::vector &w, const unsigned int n, const int a,
          const int b);

fk::matrix legendre(fk::vector const x, const unsigned int k);

fk::matrix dlegendre(fk::vector const x, const unsigned int k);

#endif
