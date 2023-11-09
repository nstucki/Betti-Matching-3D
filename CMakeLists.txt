# CMakeLists.txt

cmake_minimum_required(VERSION 3.5)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(Python_VIRTUALENV FIRST)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

project(BettiMatching)

find_package(pybind11 REQUIRED)
find_package(Threads REQUIRED)

pybind11_add_module(betti_matching
    src/_BettiMatching.cpp
    src/BettiMatching.cpp
    src/utils.cpp
    src/data_structures.cpp
    src/src_1D/BettiMatching.cpp
    src/src_1D/dimension_0.cpp
    src/src_1D/data_structures.cpp
    src/src_2D/BettiMatching.cpp
    src/src_2D/dimension_0.cpp
    src/src_2D/dimension_1.cpp
    src/src_2D/data_structures.cpp
    src/src_3D/BettiMatching.cpp
    src/src_3D/dimension_0.cpp
    src/src_3D/dimension_1.cpp
    src/src_3D/dimension_2.cpp
    src/src_3D/data_structures.cpp
    src/src_3D/enumerators.cpp
    src/src_nD/BettiMatching.cpp
    src/src_nD/dimension_0.cpp
    src/src_nD/inter_dimensions.cpp
    src/src_nD/top_dimension.cpp
    src/src_nD/data_structures.cpp
    src/src_nD/enumerators.cpp
)

# add_library(betti_matching MODULE src/BettiMatching.cpp)
# target_link_libraries(betti_matching PRIVATE pybind11::module)

add_executable(BettiMatching
    src/main.cpp
    src/BettiMatching.cpp
    src/utils.cpp
    src/data_structures.cpp
    src/src_1D/BettiMatching.cpp
    src/src_1D/dimension_0.cpp
    src/src_1D/data_structures.cpp
    src/src_2D/BettiMatching.cpp
    src/src_2D/dimension_0.cpp
    src/src_2D/dimension_1.cpp
    src/src_2D/data_structures.cpp
    src/src_3D/BettiMatching.cpp
    src/src_3D/dimension_0.cpp
    src/src_3D/dimension_1.cpp
    src/src_3D/dimension_2.cpp
    src/src_3D/data_structures.cpp
    src/src_3D/enumerators.cpp
    src/src_nD/BettiMatching.cpp
    src/src_nD/dimension_0.cpp
    src/src_nD/inter_dimensions.cpp
    src/src_nD/top_dimension.cpp
    src/src_nD/data_structures.cpp
    src/src_nD/enumerators.cpp
)

set_target_properties(BettiMatching betti_matching
PROPERTIES
CXX_STANDARD 17
CXX_EXTENSIONS OFF
)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)

target_link_libraries(BettiMatching Threads::Threads)
