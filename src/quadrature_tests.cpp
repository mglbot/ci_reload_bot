
#include "quadrature.hpp"

#include "matlab_utilities.hpp"
#include "tests_general.hpp"

TEST_CASE("lgwt() matches the matlab implementation", "[setup]")
{
  SECTION("lgwt(10, -1, 1)")
  {
    fk::vector xtest;
    fk::vector wtest;
    fk::vector xgold =
        readVectorFromTxtFile("../testing/generated-inputs/lgwt_x10.dat");
    fk::vector wgold =
        readVectorFromTxtFile("../testing/generated-inputs/lgwt_w10.dat");

    lgwt(xtest, wtest, 10, -1, 1);

    REQUIRE(xgold.size() == 10);
    REQUIRE(wgold.size() == 10);
    REQUIRE(xtest == xgold);
    REQUIRE(wtest == wgold);
  }

  SECTION("lgwt(100, -13, 34)")
  {
    fk::vector xtest;
    fk::vector wtest;
    fk::vector xgold =
        readVectorFromTxtFile("../testing/generated-inputs/lgwt_x100.dat");
    fk::vector wgold =
        readVectorFromTxtFile("../testing/generated-inputs/lgwt_w100.dat");

    lgwt(xtest, wtest, 100, -13, 34);

    REQUIRE(xgold.size() == 100);
    REQUIRE(wgold.size() == 100);
    REQUIRE(xtest == xgold);
    REQUIRE(wtest == wgold);
  }
}

TEST_CASE("[d]legendre() matches the matlab implementation on -1,1", "[setup]")
{
  fk::vector xin;
  fk::vector win;
  lgwt(xin, win, 10, -1, 1);

  SECTION("legendre(x10, 2)")
  {
    fk::matrix gold = fk::readMatrixFromTxtFile(
        "../testing/generated-inputs/legendre_10_2.dat");
    fk::matrix test = legendre(xin, 2);
    REQUIRE(test == gold.transpose());
  }

  SECTION("dlegendre(x10, 2)")
  {
    fk::matrix gold = fk::readMatrixFromTxtFile(
        "../testing/generated-inputs/dlegendre_10_2.dat");
    fk::matrix test = dlegendre(xin, 2);
    REQUIRE(test == gold.transpose());
  }
}

TEST_CASE("[d]legendre() matches the matlab implementation on -13,34",
          "[setup]")
{
  fk::vector xin;
  fk::vector win;
  lgwt(xin, win, 100, -13, 34);

  SECTION("legendre(x100, 2)")
  {
    fk::matrix gold = fk::readMatrixFromTxtFile(
        "../testing/generated-inputs/legendre_100_2.dat");
    fk::matrix test = legendre(xin, 2);
    REQUIRE(test == gold.transpose());
  }

  SECTION("dlegendre(x100, 2)")
  {
    fk::matrix gold = fk::readMatrixFromTxtFile(
        "../testing/generated-inputs/dlegendre_100_2.dat");
    fk::matrix test = dlegendre(xin, 2);
    REQUIRE(test == gold.transpose());
  }
}
