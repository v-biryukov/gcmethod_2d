project(mesh_2d)

set (${PROJECT_NAME}_headers mesh_2d.h)
set (${PROJECT_NAME}_sources mesh_2d.cpp)

include_directories (${CMAKE_SOURCE_DIR}/lib/triangle)
link_directories (${CMAKE_SOURCE_DIR}/lib/triangle)

add_library (${PROJECT_NAME} STATIC ${${PROJECT_NAME}_headers} ${${PROJECT_NAME}_sources} ${CMAKE_SOURCE_DIR}/lib/triangle/triangle.o)

