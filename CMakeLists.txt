
cmake_minimum_required(VERSION 2.8)

PROJECT(FDTD_3d)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(OpenMP)
if (OPENMP_FOUND) 
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_executable(FDTD_3d MACOSX_BUNDLE fdtd_3d.cpp init_pml.cpp)

target_link_libraries(FDTD_3d ${VTK_LIBRARIES})
