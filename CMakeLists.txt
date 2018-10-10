cmake_minimum_required (VERSION 3.9)
project (gfk6d)

###############################################################################
## Set up the compiler and general global build options
###############################################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options("-g" "-Wall")
message("CMAKE_BUILD_TYPE is ${CMAKE_BUILD_TYPE}")

# set up possible commandline input variable defaults
set(BUILD_OPENBLAS FALSE CACHE BOOL "Download and build our own OpenBLAS")
set(TRAVIS FALSE CACHE BOOL "Enable fixes for building in Travis CI VMs")

###############################################################################
## External support
## FIXME use URL_HASH to skip download
###############################################################################

# if no request to specifically build OpenBLAS, see what we can find
if( NOT BUILD_OPENBLAS )
  find_package(LAPACK)
endif()
if( LAPACK_FOUND )
  message( STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}" )
  message( STATUS "LAPACK linker flags: ${LAPACK_LINKER_FLAGS}" )
  set(LINALG_LIBS ${LAPACK_LIBRARIES})
endif()

# if cmake couldn't find other blas/lapack, or the user asked to build openblas
if( NOT LAPACK_FOUND )

  # first check if it has already been built
  set(OpenBLAS_PATH ${CMAKE_SOURCE_DIR}/contrib/blas/openblas)
  find_package(OpenBLAS PATHS ${OpenBLAS_PATH})
  if( OpenBLAS_FOUND )
    set(LINALG_LIBS ${OpenBLAS_PATH}/lib/libopenblas.so)
    message(STATUS "OpenBLAS library: ${OpenBLAS_LIBRARIES}")

  # build it if necessary
  else( NOT DEFINED OpenBLAS_LIBRARIES )
    message(STATUS "OpenBLAS not found. Will build from source.")

    include(ExternalProject)
    ExternalProject_Add(openblas-ext
      PREFIX contrib/blas/openblas
      URL https://github.com/xianyi/OpenBLAS/archive/v0.3.0.tar.gz
      CONFIGURE_COMMAND ""
      BUILD_COMMAND make
      BUILD_IN_SOURCE 1
      INSTALL_COMMAND make PREFIX=${OpenBLAS_PATH} install
      )
    set(BUILD_OpenBLAS TRUE)
  endif()
endif()

###############################################################################
## Building gfk6d
#
# to add a component:
# 1) use add_library() to define the component and its source
# 2) use target_link_libraries() to define its external component dependencies
# 3) add component to ALL_GFK6D_COMPONENTS variable in alphabetical order
# 4) add the component's test (see next section)
###############################################################################

add_library(tensors src/tensors.cpp)
if( BUILD_OpenBLAS )
  add_dependencies(tensors openblas-ext)
  set(LINALG_LIBS ${OpenBLAS_PATH}/lib/libopenblas.so)
endif()
target_link_libraries(tensors ${LINALG_LIBS})

add_library(initial_conditions src/initial_conditions.cpp)
target_link_libraries(initial_conditions quadrature pde)

add_library(pde src/pde.cpp)

add_library(matlab src/matlab_utilities.cpp)
target_link_libraries(matlab tensors)

add_library(quadrature src/quadrature.cpp)
target_link_libraries(quadrature matlab)

add_library(two_scale src/operator_two_scale.cpp)
target_link_libraries(two_scale matlab tensors)

set(ALL_GFK6D_LIBS
  initial_conditions
  matlab
  quadrature
  pde
  two_scale
  tensors
  )

add_executable(gfk6d src/main.cpp)
target_link_libraries(gfk6d ${ALL_GFK6D_LIBS})

add_library(forwardMWT src/forwardMWT.cpp)



###############################################################################
## Testing gfk6d
#
# to add a component test:
# 1) add source file to TEST_SOURCES variable in alphabetical order
###############################################################################
option(BUILD_TESTS "build tests for gfk6d" ON)
if(BUILD_TESTS)
  enable_testing()

  # Prepare "Catch" library for other executables
  set(CATCH_INCLUDE_DIR contrib/catch2/include)
  add_library(Catch INTERFACE)
  target_include_directories(Catch INTERFACE
    ${CMAKE_SOURCE_DIR}/testing
    ${CATCH_INCLUDE_DIR}
    )

  # to help build within travis containers
  if( TRAVIS )
    target_compile_options(Catch INTERFACE -DCATCH_CONFIG_NO_CPP17_UNCAUGHT_EXCEPTIONS)
  endif()

  # Make test executable
  set(TEST_SOURCES
    testing/tests_general.cpp
    src/initial_conditions_tests.cpp
    src/matlab_utilities_tests.cpp
    src/operator_two_scale_tests.cpp
    src/quadrature_tests.cpp
    src/pde.cpp
    src/tensors_tests.cpp
    src/forwardMWT_tests.cpp
    )

  add_executable(tests ${TEST_SOURCES})
  target_link_libraries(tests ${ALL_GFK6D_LIBS} Catch)

  add_test(all-tests tests)
endif()
