cmake_minimum_required(VERSION 3.26)
project(FluidSim)

set(CMAKE_CXX_STANDARD 17)

# Set paths to your libraries
set(GLFW_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glfw-3.4)
set(GLEW_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glew-2.1.0)
set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/imgui)

# Build GLFW from source
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory(${GLFW_DIR})

# Include directories
include_directories(
        ${GLFW_DIR}/include
        ${GLEW_DIR}/include
        ${IMGUI_DIR}
        ${IMGUI_DIR}/backends
)

# ImGui source files
set(IMGUI_SOURCES
        ${IMGUI_DIR}/imgui.cpp
        ${IMGUI_DIR}/imgui_draw.cpp
        ${IMGUI_DIR}/imgui_tables.cpp
        ${IMGUI_DIR}/imgui_widgets.cpp
        ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
        ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
)

# GLEW source files (build from source)
set(GLEW_SOURCES
        ${GLEW_DIR}/src/glew.c
)

# Define GLEW_STATIC to use static linking
add_definitions(-DGLEW_STATIC)

# Your application
add_executable(FluidSim
        main.cpp
        ${IMGUI_SOURCES}
        ${GLEW_SOURCES}
        Grid.h
        Grid.cpp
)

# Link libraries
target_link_libraries(FluidSim
        glfw
        opengl32
        gdi32
)