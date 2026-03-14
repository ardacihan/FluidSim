cmake_minimum_required(VERSION 3.26)
project(FluidSim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Paths
set(GLFW_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glfw-3.4)
set(GLEW_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glew-2.1.0)
set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/imgui)

# Build GLFW
if(EXISTS ${GLFW_DIR}/CMakeLists.txt)
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    add_subdirectory(${GLFW_DIR})
    set(GLFW_LIB glfw)
else()
    message(FATAL_ERROR "GLFW not found")
endif()

# Include dirs
include_directories(
        ${GLFW_DIR}/include
        ${GLEW_DIR}/include
        ${IMGUI_DIR}
        ${IMGUI_DIR}/backends
        ${CMAKE_CURRENT_SOURCE_DIR}
)

# ImGui sources
set(IMGUI_SOURCES
        ${IMGUI_DIR}/imgui.cpp
        ${IMGUI_DIR}/imgui_draw.cpp
        ${IMGUI_DIR}/imgui_tables.cpp
        ${IMGUI_DIR}/imgui_widgets.cpp
        ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
        ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
)

# GLEW
set(GLEW_SOURCES ${GLEW_DIR}/src/glew.c)
add_definitions(-DGLEW_STATIC)

# Executable
add_executable(FluidSim
        main.cpp
        main2d.cpp
        main2d.h
        Grid2D.h
        ${IMGUI_SOURCES}
        ${GLEW_SOURCES}
)

# Link libraries
target_link_libraries(FluidSim
        ${GLFW_LIB}
        opengl32
        glu32
        gdi32
        user32
        shell32
)

# ImGui OpenGL backend
target_compile_definitions(FluidSim PRIVATE
        IMGUI_IMPL_OPENGL_LOADER_GLEW
        GLEW_STATIC
)
