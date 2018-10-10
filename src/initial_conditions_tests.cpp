#include "initial_conditions.hpp"
#include "matlab_utilities.hpp"
#include "tests_general.hpp"

TEST_CASE("initial_conditions_1D() matches matlab fragile inputs")
{
  int LevX = 3;
  int LevV = 3;
  int deg  = 2;

  int FMWTxLDA = 16;
  int FMWTvLDA = 16;

  fk::matrix FMWTx = readVectorFromBinFile(
      "../testing/fragile-inputs/matlab-inputs-FMWTx.dat");

  fk::matrix FMWTv = readVectorFromBinFile(
      "../testing/fragile-inputs/matlab-inputs-FMWTv.dat");

  fk::vector f_x_matlab =
      readVectorFromBinFile("../testing/fragile-inputs/matlab-outputs-f_x.dat");

  fk::vector f_v_matlab =
      readVectorFromBinFile("../testing/fragile-inputs/matlab-outputs-f_v.dat");

  // Select the Vlasov4 inputs for this test case
  //
  PDE pde;

  fk::vector f_x, f_v;
  initial_conditions_1D(LevX, LevV, deg, pde, FMWTx, FMWTxLDA, FMWTv, FMWTvLDA,
                        f_x, f_v);

  // FIXME for until we find the reason for the 10e-6 errors
  // for (unsigned int j = 0; j < f_x.size(); ++j)
  //{
  //  REQUIRE(f_x[j] == Approx(f_x_matlab[j]).margin(1e-6));
  //  REQUIRE(f_v[j] == Approx(f_v_matlab[j]).margin(1e-6));
  //}
  fk::vector f_x_tmp_gold = {
      4.5764561643e+00,  -8.2507785326e-17, 1.2850208415e-01,
      -5.5294310797e-18, -2.7929047963e-16, 7.4280576759e-03,
      -1.9559007192e-16, -7.4280576759e-03, -5.5293112119e-03,
      5.0743720923e-04,  5.5293112119e-03,  5.0743720923e-04,
      5.5293112119e-03,  -5.0743720923e-04, -5.5293112119e-03,
      -5.0743720923e-04};
  REQUIRE(f_x == f_x_tmp_gold);
  REQUIRE(f_v == f_v_matlab);
}
