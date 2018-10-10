#include "matlab_utilities.hpp"
#include "operator_two_scale.hpp"
#include "tests_general.hpp"
#include <iostream>
#include <vector>

TEST_CASE("operator_two_scale function working appropriately", "[operators]")
{
  SECTION("operator_two_scale(2, 3)")
  {
    fk::matrix gold = fk::readMatrixFromTxtFile(
        "../testing/generated-inputs/operator_two_scale_2_3.dat");
    fk::matrix test = operator_two_scale(2, 3);
    REQUIRE(gold == test);
  }
  SECTION("operator_two_scale(3, 4)")
  {
    fk::matrix gold = fk::readMatrixFromTxtFile(
        "../testing/generated-inputs/operator_two_scale_3_4.dat");
    fk::matrix test = operator_two_scale(3, 4);
    REQUIRE(gold == test);
  }
  SECTION("operator_two_scale(4, 5)")
  {
    fk::matrix gold = fk::readMatrixFromTxtFile(
        "../testing/generated-inputs/operator_two_scale_4_5.dat");
    fk::matrix test = operator_two_scale(4, 5);
    REQUIRE(gold == test);
  }
  SECTION("operator_two_scale(2, 6)")
  {
    fk::matrix gold = fk::readMatrixFromTxtFile(
        "../testing/generated-inputs/operator_two_scale_2_6.dat");
    fk::matrix test = operator_two_scale(2, 6);
    REQUIRE(gold == test);
  }
}
