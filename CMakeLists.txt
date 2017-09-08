
cmake_minimum_required(VERSION 2.8)

PROJECT(FDTD_3d)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(FDTD_3d MACOSX_BUNDLE fdtd_3d.cpp )

target_link_libraries(FDTD_3d ${VTK_LIBRARIES})