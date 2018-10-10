#include "operator_two_scale.hpp"
#include "matlab_utilities.hpp"
#include "tensors.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

fk::matrix operator_two_scale(int const maxDeg, int const maxLev)
{
  // Load .dat files, clean up the data, and transform. store in const vars
  auto const [g0, g1, h0, h1] = [maxDeg]() -> std::array<fk::vector, 4> {
    // FIXME it would be nicer to just generate this data on the fly
    fk::vector g0(readVectorFromBinFile("../data/two_scale_rel_G0_" +
                                        std::to_string(maxDeg) + ".dat"));
    fk::vector h0(readVectorFromBinFile("../data/two_scale_rel_H0_" +
                                        std::to_string(maxDeg) + ".dat"));

    for (auto i = 0; i < g0.size(); ++i)
    {
      if (std::abs(g0(i)) < 1e-5) g0(i) = 0.0; // clean up the data
      if (std::abs(h0(i)) < 1e-5) h0(i) = 0.0;
    }

    // now g1, h1 from g0, h0
    fk::vector g1(maxDeg * maxDeg);
    fk::vector h1(maxDeg * maxDeg);
    for (int j_x = 1; j_x <= maxDeg; j_x++)
      for (int j_y = 1; j_y <= maxDeg; j_y++)
      {
        auto idx = (j_x - 1.0) * maxDeg + (j_y - 1.0);
        g1(idx)  = std::pow(-1.0, (maxDeg + j_x + j_y - 2.0)) * g0(idx);
        h1(idx)  = std::pow(-1.0, (j_x + j_y - 2.0)) * h0(idx);
      }

    return {g0, g1, h0, h1};
  }();

  fk::matrix fmwt(maxDeg * maxLev, maxDeg * maxLev);

  for (int j = 1; j <= (float)maxLev / 2.0; j++)
  {
    int rs = maxDeg * (j - 1.0) + 1.0;
    int cs = 2.0 * maxDeg * (j - 1.0) + 1.0;

    for (int j_x = 1; j_x <= maxDeg; j_x++)
    {
      for (int j_y = 1; j_y <= maxDeg; j_y++)
      {
        int const row1 = rs + j_x - 1 - 1;
        int col1       = cs + j_y - 1 - 1;
        int const row2 = j_x - 1;
        int const col2 = j_y - 1;

        fmwt(row1, col1) = h0(col2 * maxDeg + row2);

        col1             = cs + maxDeg + j_y - 1 - 1;
        fmwt(row1, col1) = h1(col2 * maxDeg + row2);
      }
    }

    rs = maxDeg * (j + (float)maxLev / 2.0 - 1.0) + 1.0;
    cs = 2.0 * maxDeg * (j - 1.0) + 1.0;

    for (int j_x = 1; j_x <= maxDeg; j_x++)
    {
      for (int j_y = 1; j_y <= maxDeg; j_y++)
      {
        int const row1 = rs + j_x - 1 - 1;
        int col1       = cs + j_y - 1 - 1;
        int const row2 = j_x - 1;
        int const col2 = j_y - 1;

        fmwt(row1, col1) = g0(col2 * maxDeg + row2);

        col1             = cs + maxDeg + j_y - 1 - 1;
        fmwt(row1, col1) = g1(col2 * maxDeg + row2);
      }
    }
  }

  fk::matrix fmwt_comp = eye(maxDeg * maxLev, maxDeg * maxLev);
  fk::matrix cfmwt(maxDeg * maxLev, maxDeg * maxLev);

  int const n = std::floor(std::log2(maxLev));
  for (int j = 1; j <= n; j++)
  {
    cfmwt = fmwt;

    if (j > 1)
    {
      cfmwt = fk::matrix(cfmwt.nrows(), cfmwt.ncols()); // zero out

      int const cn = std::pow(2.0, n - j + 1.0) * maxDeg;

      int rs       = cn + 1;
      int const cs = cn + 1;

      for (int ii = 0; ii <= maxDeg * maxLev - (cn + 1); ii++)
        for (int jj = 0; jj <= maxDeg * maxLev - (cn + 1); jj++)
          if (ii == jj)
          {
            int const row = rs + ii - 1;
            int const col = cs + jj - 1;

            cfmwt(row, col) = 1;
          }

      for (int ii = 1; ii <= (float)cn / 2.0; ii++)
        for (int jj = 1; jj <= cn; jj++)
        {
          int const row = ii - 1;
          int const col = jj - 1;

          cfmwt(row, col) = fmwt(row, col);
        }

      rs = (float)(maxDeg * maxLev) / 2.0 + 1;
      for (int ii = 0; ii <= (float)cn / 2.0 - 1; ii++)
        for (int jj = 1; jj <= cn; jj++)
        {
          int const row1 = (float)cn / 2.0 + 1 + ii - 1;
          int const col1 = jj - 1;

          int const row2 = rs + ii - 1;
          int const col2 = jj - 1;

          cfmwt(row1, col1) = fmwt(row2, col2);
        }
    }

    fk::matrix tmp = cfmwt * fmwt_comp;
    fmwt_comp      = tmp;
  }

  return fmwt_comp;
}
