''FP Test'' is a small C++11 program which tests elementary math
functions in C standard library implementation.

Requirements
------------

- C++11 compatible compiler: [g++](http://gcc.gnu.org/) (version >= 4.8.1), or [clang++](http://clang.llvm.org/cxx_status.html) (version >= 3.3)
- [CMake](http://www.cmake.org)
- [GMP (GNU multiprecision library)](http://gmplib.org/)
- [MPFR (GNU MPFR Library)](http://www.mpfr.org/)

Build
-----
Instructions for DEBUG build

    mkdir -p build/debug
    cd build/debug
    cmake -DCMAKE_CXX_COMPILER=<c++-compiler> -DCMAKE_BUILD_TYPE=DEBUG ../../src
    make

Instructions for RELEASE build

    mkdir -p build/release
    cd build/release
    cmake -DCMAKE_CXX_COMPILER=<c++-compiler> -DCMAKE_BUILD_TYPE=RELEASE ../../src
    make
