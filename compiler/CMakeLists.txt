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
		glew
 )

if (APPLE)
    find_package (OpenGL REQUIRED)
    include_directories(/usr/local/include)
    find_library(COCOA_LIBRARY Cocoa REQUIRED)
    find_library(IOKIT_LIBRARY IOKit REQUIRED)
    find_library(COREVID_LIBRARY CoreVideo REQUIRED)
endif (APPLE)

set (HEADERS
		../src/Solver.h
		)
set (SOURCES 
		../src/Solver.cpp
		../src/main.cpp)

link_directories (
		../3rdparty/GLEW/lib/
		../3rdparty/GLFW/lib/
		../3rdparty/boost/lib
		)

source_group ("Header Files" FILES ${HEADERS})
source_group ("Source Files" FILES ${SOURCES})

set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

add_executable (${PROJECT} ${HEADERS} ${SOURCES})
target_link_libraries (${PROJECT} ${OPENGL_LIBRARIES} ${LIBRARIES} ${COCOA_LIBRARY} ${COREVID_LIBRARY} ${IOKIT_LIBRARY})
