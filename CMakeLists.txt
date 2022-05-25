cmake_minimum_required(VERSION 3.2.0)
project(bsm VERSION 0.1.0)
set (CMAKE_CXX_STANDARD 17)

add_library(gnuplotcpp INTERFACE)
add_library(tester STATIC incl/test.hpp src/test.cpp)
add_library(psr STATIC incl/psr.hpp src/psr.cpp)
add_library(kdt STATIC incl/kdt.hpp src/kdt.cpp)

add_library(asm STATIC incl/asm.hpp src/asm.cpp)


include(CTest)
enable_testing()

add_executable(bsm main.cpp)

target_link_libraries(bsm PRIVATE tester)
target_link_libraries(bsm PRIVATE psr)
target_link_libraries(bsm PRIVATE kdt)

target_link_libraries(bsm PRIVATE asm)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
