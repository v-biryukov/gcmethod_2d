cmake_minimum_required (VERSION 2.6)
project (gcmethod_2d)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(${PROJECT_NAME}_headers src/convection_equation_solver.h src/linear_elastisity_2d.h src/tensor2d.h src/vector2d.h)
set(${PROJECT_NAME}_sources src/convection_equation_solver.cpp src/linear_elastisity_2d.cpp src/main.cpp)


add_library (triangle ${CMAKE_SOURCE_DIR}/lib/triangle/triangle.c) 
include_directories (${CMAKE_SOURCE_DIR}/lib/triangle)
link_directories (${CMAKE_SOURCE_DIR}/lib/triangle)

add_executable (${PROJECT_NAME} ${${PROJECT_NAME}_headers} ${${PROJECT_NAME}_sources})
target_link_libraries (${PROJECT_NAME} triangle mesh_2d)

add_subdirectory (src/mesh_2d)
