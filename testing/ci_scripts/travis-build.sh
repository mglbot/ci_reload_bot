#!/bin/bash

echo ""
echo ""
echo "starting new test for real full precision"
echo ""
echo ""

STATUS=0

env

$CC --version
$CXX --version

mkdir travis_build
cd travis_build

cmake .. -DTRAVIS=1

make

echo
echo checking all tests
echo ----------------------------------------------------
echo

./tests
