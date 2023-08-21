# CMakeLists.txt

cmake_minimum_required(VERSION 3.5)
project(bettiMatching)

#find_package(pybind11 REQUIRED)
#pybind11_add_module(betti_matching
#    src/src_nD/_betti_matching_nD.cpp
#    src/src_nD/betti_matching_nD.cpp
#    src/src_nD/BettiMatching.cpp
#    src/src_nD/dimension_0.cpp
#    src/src_nD/inter_dimensions.cpp
#    src/src_nD/top_dimension.cpp
#    src/src_nD/data_structures.cpp
#    src/src_nD/enumerators.cpp
#    src/src_nD/utils.cpp
#)

add_executable(BettiMatching
    src/main.cpp
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

set_target_properties(BettiMatching
PROPERTIES
CXX_STANDARD 17
CXX_EXTENSIONS OFF
)
