cmake_minimum_required(VERSION 3.15)
project(de_c_pure)

find_package(Armadillo REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 17)

add_executable(de_c_pure main.cpp lattice_base.h analysis.h)
target_link_libraries(de_c_pure ${LAPACK_LIBRARIES} ${BLAS_LIBARIES} ${ARMADILLO_LIBRARIES} ${TORCH_LIBRARIES})