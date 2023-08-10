# CMakeLists.txt
set(CMAKE_CXX_FLAGS
"-include /usr/include/c++/13.1.1/limits \
-include /usr/include/c++/13.1.1/cassert \
-include /usr/include/c++/13.1.1/algorithm \
-include /usr/include/c++/13.1.1/cstdint \
-include /usr/include/c++/13.1.1/utility \
-I usr/lib/python3.11/site-packages/pybind11/include")

set(CMAKE_BUILD_TYPE, "Release")

cmake_minimum_required(VERSION 2.8.12)
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

add_executable(BettiMatching
    src/src_nD/BettiMatching.cpp
    src/src_nD/dimension_0.cpp
    src/src_nD/inter_dimensions.cpp
    src/src_nD/top_dimension.cpp
    src/src_nD/data_structures.cpp
    src/src_nD/enumerators.cpp
    src/src_nD/utils.cpp
    src/src_nD/betti_matching_nD.cpp
)

# set(CPLUS_INCLUDE_PATH "/usr/include/c++/13.1.1")

# target_compile_options(BettiMatching "-include /usr/include/c++/13.1.1/limits")
#  -include /usr/include/c++/13.1.1/cassert -include /usr/include/c++/13.1.1/algorithm -include /usr/include/c++/13.1.1/cstdint)