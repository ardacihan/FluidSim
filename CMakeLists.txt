cmake_minimum_required(VERSION 3.26)
project(FluidSim)

set(CMAKE_CXX_STANDARD 17)

# Set paths to your libraries
set(GLFW_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glfw-3.4)
set(GLEW_DIR ${CMAKE_CURRENT_SOURCE_DIR}/glew-2.1.0)
set(IMGUI_DIR ${CMAKE_CURRENT_SOURCE_DIR}/imgui)

# First, check if GLFW has its own CMakeLists.txt
if(EXISTS ${GLFW_DIR}/CMakeLists.txt)
    # Build GLFW as a subdirectory project
    set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    add_subdirectory(${GLFW_DIR})

    # GLFW target will be named 'glfw' or 'glfw3'
    set(GLFW_LIB glfw)
else()
    message(FATAL_ERROR "GLFW not found or doesn't have CMakeLists.txt. Please download GLFW from https://www.glfw.org/download.html")
endif()

# Include directories
include_directories(
        ${GLFW_DIR}/include
        ${GLEW_DIR}/include
        ${IMGUI_DIR}
        ${IMGUI_DIR}/backends
        ${CMAKE_CURRENT_SOURCE_DIR} # Include current dir for Grid3D.h
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

# GLEW sources - compile it directly
set(GLEW_SOURCES ${GLEW_DIR}/src/glew.c)
add_definitions(-DGLEW_STATIC)

# Executable
add_executable(FluidSim
        main.cpp
        ${IMGUI_SOURCES}
        ${GLEW_SOURCES}
        Grid3D.cpp
        main2d.cpp
        main2d.h
)

# Link libraries
target_link_libraries(FluidSim
        ${GLFW_LIB}        # GLFW library
        opengl32          # Windows OpenGL
        glu32             # GLU utilities (optional)
)

# Additional Windows libraries needed by GLFW
target_link_libraries(FluidSim
        gdi32
        user32
        shell32
)

# For ImGui OpenGL 3 backend
target_compile_definitions(FluidSim PRIVATE
        IMGUI_IMPL_OPENGL_LOADER_GLEW
        GLEW_STATIC
)

# Optional: Set properties for the target
set_target_properties(FluidSim PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
)