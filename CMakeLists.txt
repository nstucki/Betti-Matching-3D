# CMakeLists.txt

set(CMAKE_BUILD_TYPE, "Release")

cmake_minimum_required(VERSION 3.5)
project(bettiMatching)

find_package(pybind11 REQUIRED)
pybind11_add_module(betti_matching
    src/src_nD/_betti_matching_nD.cpp
    src/src_nD/betti_matching_nD.cpp
    src/src_nD/BettiMatching.cpp
    src/src_nD/dimension_0.cpp
    src/src_nD/inter_dimensions.cpp
    src/src_nD/top_dimension.cpp
    src/src_nD/data_structures.cpp
    src/src_nD/enumerators.cpp
    src/src_nD/utils.cpp
)

add_executable(BettiMatchingND
    src/src_nD/BettiMatching.cpp
    src/src_nD/dimension_0.cpp
    src/src_nD/inter_dimensions.cpp
    src/src_nD/top_dimension.cpp
    src/src_nD/data_structures.cpp
    src/src_nD/enumerators.cpp
    src/src_nD/utils.cpp
    src/src_nD/betti_matching_nD.cpp
)

add_executable(BettiMatching
    src/BettiMatching.cpp
    src/utils.cpp
    src/src_3D/BettiMatching3D.cpp
    src/src_3D/dimension_0.cpp
    src/src_3D/dimension_1.cpp
    src/src_3D/dimension_2.cpp
    src/src_3D/data_structures.cpp
    src/src_3D/enumerators.cpp
)

set_target_properties(BettiMatchingND BettiMatching betti_matching
PROPERTIES
CXX_STANDARD 17
CXX_EXTENSIONS OFF
)
