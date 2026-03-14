# FluidSim

Real-time 2D and 3D fluid simulation using a staggered MAC grid, semi-Lagrangian advection, and pressure projection.

## Dependencies

- [GLFW 3.4](https://www.glfw.org/download.html)
- [GLEW 2.1.0](https://sourceforge.net/projects/glew/)
- [ImGui](https://github.com/ocornut/imgui) (master branch)

Place each library in the project root alongside `CMakeLists.txt`.

## Switching Between 2D and 3D

**2D** (`main2d.cpp`) — `main()` is defined here. Runs by default.

**3D** (`main.cpp`) — entry point is named `main2()`. To run the 3D version, rename it to `main()` and rename the 2D `main()` to something else (e.g. `main2d()`).

## Controls

| Action | Input |
|---|---|
| Add velocity / dye | `Ctrl + Left Click` |
| Pan (2D) | Middle Mouse |
| Zoom (2D) | Scroll Wheel |
| Orbit (3D) | Right Mouse Drag |
