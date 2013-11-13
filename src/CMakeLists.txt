cmake_minimum_required(VERSION 2.8.7)
project(FP_BUGS_TEST CXX)
set(CMAKE_COLOR_MAKEFILE ON)

# Set Module Path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules/")

# Initialize CXXFLAGS.
set(CMAKE_CXX_FLAGS                "-Wall -Wextra -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG          "-g -DDEBUG -DTRACE")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

# MPFR
find_package(MPFR 3.1.0)
include_directories(${MPFR_INCLUDES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${MPFR_LIBRARIES})

# GMP
find_package(GMP 5.0.5)
include_directories(${GMP_INCLUDE_DIR})
set(EXTRA_LIBS ${EXTRA_LIBS} ${GMP_LIBRARIES})

add_executable(fp_bugs_test fp_bugs_test.cpp)
target_link_libraries(fp_bugs_test ${EXTRA_LIBS} -static-libgcc -static-libstdc++)
