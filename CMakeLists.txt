# CMakeLists.txt

cmake_minimum_required(VERSION 3.11)
project(bettiMatching)

find_package(pybind11 REQUIRED)
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

set_target_properties(BettiMatching betti_matching
PROPERTIES
CXX_STANDARD 17
CXX_EXTENSIONS OFF
)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

set(TRACY_ENABLE 1)
set(TRACY_ON_DEMAND 0)
set(TRACY_SAMPLING_HZ 2000)
set(TRACY_NO_EXIT 0)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
if(COMPILER_SUPPORTS_MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

include(FetchContent)
FetchContent_Declare (
    tracy
    GIT_REPOSITORY https://github.com/wolfpld/tracy.git
    GIT_TAG v0.9.1 
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable (tracy)

target_link_libraries (betti_matching PUBLIC TracyClient)
target_link_libraries (BettiMatching PUBLIC TracyClient)
set_property(TARGET betti_matching PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET BettiMatching PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET TracyClient PROPERTY POSITION_INDEPENDENT_CODE ON)
