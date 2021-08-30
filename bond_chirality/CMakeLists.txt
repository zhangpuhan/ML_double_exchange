cmake_minimum_required(VERSION 3.15)
project(de_c_torch)

# set(Torch_DIR libtorch/share/cmake/Torch)

find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
# find_package(Torch REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 20)

add_executable(de_c_torch main.cpp lattice_base.h help_function.h analysis.h real_dynamics.h)

target_link_libraries(de_c_torch ${LAPACK_LIBRARIES} ${BLAS_LIBARIES} ${ARMADILLO_LIBRARIES})