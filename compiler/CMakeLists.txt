cmake_minimum_required(VERSION 2.8)
set (PROJECT Fluid)
 
project (${PROJECT})
 
include_directories (
		../src/
		../3rdparty/boost/include/
		../3rdparty/GLEW/include/
		../3rdparty/GLFW/include/
		../3rdparty/glm/include/
		)

set (LIBRARIES
		glfw3
		glew32
		glu32
		opengl32 )

set (HEADERS
		../src/aabb.h
		../src/IRender.h
		../src/OGLRender.h
		../src/Ray.h
		../src/Solver.h
		)
set (SOURCES 
		../src/OGLRender.cpp
		../src/Solver.cpp
		../src/main.cpp)

link_directories (
		../3rdparty/GLEW/lib/
		../3rdparty/GLFW/lib/
		../3rdparty/boost/lib
		)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

add_executable (${PROJECT} ${HEADERS} ${SOURCES})
target_link_libraries (${PROJECT} ${LIBRARIES})
