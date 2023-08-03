# CMakeLists.txt
set(CMAKE_CXX_FLAGS
"-include /usr/include/c++/13.1.1/limits \
-include /usr/include/c++/13.1.1/cassert \
-include /usr/include/c++/13.1.1/algorithm \
-include /usr/include/c++/13.1.1/cstdint \
-include /usr/include/c++/13.1.1/utility \
-I usr/lib/python3.11/site-packages/pybind11/include")


cmake_minimum_required(VERSION 2.8.12)
project(bettiMatching)

add_executable(BettiMatching
    src/src_nD/BettiMatching.cpp
    src/src_nD/dimension_0.cpp
    src/src_nD/inter_dimensions.cpp
    src/src_nD/top_dimension.cpp
    src/src_nD/data_structures.cpp
    src/src_nD/enumerators.cpp
    src/src_nD/utils.cpp
)

# set(CPLUS_INCLUDE_PATH "/usr/include/c++/13.1.1")

# target_compile_options(BettiMatching "-include /usr/include/c++/13.1.1/limits")
#  -include /usr/include/c++/13.1.1/cassert -include /usr/include/c++/13.1.1/algorithm -include /usr/include/c++/13.1.1/cstdint)