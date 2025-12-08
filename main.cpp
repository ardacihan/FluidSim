#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "Grid3D.h"

// --- Configuration ---
struct Config {
    static const int GRID_SIZE = 10;
    static const int WINDOW_WIDTH = 1280;
    static const int WINDOW_HEIGHT = 720;
};

// --- Global State ---
struct CameraState {
    float rotX = 30.0f;
    float rotY = -45.0f;
    float distance = 80.0f;
    ImVec2 lastMousePos;
    bool isDragging = false;
};

struct SimulationState {
    float dt = 0.1f;
    float diff = 0.0001f;
    float visc = 0.0001f;
    bool isRunning = true;
    bool showDivergence = false;
    bool showVelocityVectors = true;
    bool showInterpolationVectors = false;
    float vectorScale = 2.0f;
    int vectorSkip = 2; // Show every Nth vector
    float minVelocityThreshold = 0.01f;
    float interpolationOffset = 0.5f; // Offset for interpolation test points
    int interpolationSubdivisions = 5; // Number of interpolation points between grid cells
};

// --- Forward Declarations ---
namespace Graphics {
    void drawCube(float x, float y, float z, float size, float r, float g, float b, float alpha);
    void drawWireframeBox(float size);
    void drawArrow(float x1, float y1, float z1, float x2, float y2, float z2, float r, float g, float b);
    void setupPerspective(int width, int height);
    void setupCamera(const CameraState& camera);
    void renderScene(Grid3D& grid, const SimulationState& state);
}

namespace UI {
    void setupImGui(GLFWwindow* window);
    void renderImGui(Grid3D& grid, SimulationState& state, CameraState& camera);
}

namespace Input {
    void handleCameraInput(CameraState& camera);
}

namespace App {
    GLFWwindow* initializeWindow();
    void initializeGraphics();
    void mainLoop(GLFWwindow* window, Grid3D& grid, CameraState& camera, SimulationState& state);
    void shutdown(GLFWwindow* window);
}

// --- Graphics Implementation ---
void Graphics::drawCube(float x, float y, float z, float size, float r, float g, float b, float alpha) {
    glPushMatrix();
    glTranslatef(x, y, z);
    glScalef(size, size, size);
    glColor4f(r, g, b, alpha);

    glBegin(GL_QUADS);
    // Front
    glVertex3f(-0.5f, -0.5f,  0.5f); glVertex3f( 0.5f, -0.5f,  0.5f);
    glVertex3f( 0.5f,  0.5f,  0.5f); glVertex3f(-0.5f,  0.5f,  0.5f);
    // Back
    glVertex3f(-0.5f, -0.5f, -0.5f); glVertex3f(-0.5f,  0.5f, -0.5f);
    glVertex3f( 0.5f,  0.5f, -0.5f); glVertex3f( 0.5f, -0.5f, -0.5f);
    // Left
    glVertex3f(-0.5f, -0.5f, -0.5f); glVertex3f(-0.5f, -0.5f,  0.5f);
    glVertex3f(-0.5f,  0.5f,  0.5f); glVertex3f(-0.5f,  0.5f, -0.5f);
    // Right
    glVertex3f( 0.5f, -0.5f, -0.5f); glVertex3f( 0.5f,  0.5f, -0.5f);
    glVertex3f( 0.5f,  0.5f,  0.5f); glVertex3f( 0.5f, -0.5f,  0.5f);
    // Top
    glVertex3f(-0.5f,  0.5f, -0.5f); glVertex3f(-0.5f,  0.5f,  0.5f);
    glVertex3f( 0.5f,  0.5f,  0.5f); glVertex3f( 0.5f,  0.5f, -0.5f);
    // Bottom
    glVertex3f(-0.5f, -0.5f, -0.5f); glVertex3f( 0.5f, -0.5f, -0.5f);
    glVertex3f( 0.5f, -0.5f,  0.5f); glVertex3f(-0.5f, -0.5f,  0.5f);
    glEnd();

    glPopMatrix();
}

void Graphics::drawWireframeBox(float size) {
    glPushMatrix();
    glScalef(size, size, size);
    glColor3f(1.0f, 1.0f, 1.0f);
    glLineWidth(2.0f);

    glBegin(GL_LINES);
    // Bottom rect
    glVertex3f(0,0,0); glVertex3f(1,0,0);
    glVertex3f(1,0,0); glVertex3f(1,0,1);
    glVertex3f(1,0,1); glVertex3f(0,0,1);
    glVertex3f(0,0,1); glVertex3f(0,0,0);
    // Top rect
    glVertex3f(0,1,0); glVertex3f(1,1,0);
    glVertex3f(1,1,0); glVertex3f(1,1,1);
    glVertex3f(1,1,1); glVertex3f(0,1,1);
    glVertex3f(0,1,1); glVertex3f(0,1,0);
    // Vertical pillars
    glVertex3f(0,0,0); glVertex3f(0,1,0);
    glVertex3f(1,0,0); glVertex3f(1,1,0);
    glVertex3f(1,0,1); glVertex3f(1,1,1);
    glVertex3f(0,0,1); glVertex3f(0,1,1);
    glEnd();

    glPopMatrix();
}

void Graphics::drawArrow(float x1, float y1, float z1, float x2, float y2, float z2, float r, float g, float b) {
    // Calculate direction
    float dx = x2 - x1;
    float dy = y2 - y1;
    float dz = z2 - z1;
    float length = std::sqrt(dx*dx + dy*dy + dz*dz);

    if (length < 0.001f) return; // Too short to draw

    // Normalize
    dx /= length;
    dy /= length;
    dz /= length;

    // Arrow parameters
    float arrowHeadLength = std::min(length * 0.3f, 0.5f);
    float arrowHeadWidth = arrowHeadLength * 0.4f;

    // Calculate arrow head base position
    float baseX = x2 - dx * arrowHeadLength;
    float baseY = y2 - dy * arrowHeadLength;
    float baseZ = z2 - dz * arrowHeadLength;

    glDisable(GL_LIGHTING);
    glColor3f(r, g, b);
    glLineWidth(2.0f);

    // Draw shaft
    glBegin(GL_LINES);
    glVertex3f(x1, y1, z1);
    glVertex3f(baseX, baseY, baseZ);
    glEnd();

    // Draw arrow head (simplified cone as lines)
    // Find perpendicular vectors for the cone base
    float perpX, perpY, perpZ;
    if (std::abs(dx) < 0.9f) {
        perpX = dy;
        perpY = -dx;
        perpZ = 0;
    } else {
        perpX = 0;
        perpY = dz;
        perpZ = -dy;
    }
    float perpLen = std::sqrt(perpX*perpX + perpY*perpY + perpZ*perpZ);
    perpX = (perpX / perpLen) * arrowHeadWidth;
    perpY = (perpY / perpLen) * arrowHeadWidth;
    perpZ = (perpZ / perpLen) * arrowHeadWidth;

    // Draw cone approximation (4 triangular faces)
    glBegin(GL_LINES);
    for (int i = 0; i < 4; i++) {
        float angle = i * 3.14159f / 2.0f;
        float cos_a = std::cos(angle);
        float sin_a = std::sin(angle);

        // Rotate perpendicular vector
        float px = baseX + perpX * cos_a;
        float py = baseY + perpY * cos_a;
        float pz = baseZ + perpZ * sin_a;

        // Line from base point to tip
        glVertex3f(px, py, pz);
        glVertex3f(x2, y2, z2);

        // Line to next base point
        float angle2 = (i + 1) * 3.14159f / 2.0f;
        float cos_a2 = std::cos(angle2);
        float sin_a2 = std::sin(angle2);
        float px2 = baseX + perpX * cos_a2;
        float py2 = baseY + perpY * cos_a2;
        float pz2 = baseZ + perpZ * sin_a2;

        glVertex3f(px, py, pz);
        glVertex3f(px2, py2, pz2);
    }
    glEnd();

    glEnable(GL_LIGHTING);
}

void Graphics::setupPerspective(int width, int height) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    float aspect = (float)width / (float)height;
    float nearPlane = 1.0f;
    float farPlane = 500.0f;
    float fieldOfView = 60.0f;
    float fH = tan(fieldOfView / 360.0f * 3.14159f) * nearPlane;
    float fW = fH * aspect;

    glFrustum(-fW, fW, -fH, fH, nearPlane, farPlane);
}

void Graphics::setupCamera(const CameraState& camera) {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glTranslatef(0.0f, 0.0f, -camera.distance);
    glRotatef(camera.rotX, 1.0f, 0.0f, 0.0f);
    glRotatef(camera.rotY, 0.0f, 1.0f, 0.0f);
    glTranslatef(-(float)Config::GRID_SIZE/2.0f, -(float)Config::GRID_SIZE/2.0f, -(float)Config::GRID_SIZE/2.0f);
}

void Graphics::renderScene(Grid3D& grid, const SimulationState& state) {
    // Draw Grid Boundary
    drawWireframeBox((float)Config::GRID_SIZE);

    // Draw Voxels
    for(int k = 1; k <= grid.N; k++) {
        for(int j = 1; j <= grid.N; j++) {
            for(int i = 1; i <= grid.N; i++) {
                if (state.showDivergence) {
                    // Draw divergence visualization
                    std::vector<float> color = grid.getDivergenceColor(i, j, k, 0.01f);
                    if (color[3] > 0.01f) {
                        drawCube((float)i + 0.5f, (float)j + 0.5f, (float)k + 0.5f,
                                0.85f, color[0], color[1], color[2], color[3]);
                    }
                } else {
                    // Draw density visualization
                    float d = grid.dens[grid.P_IX(i, j, k)];
                    if (d > 0.05f) {
                        float alpha = std::min(d/100.0f, 0.8f);
                        float blue = 1.0f;
                        float green = std::min(d/100.0f, 1.0f);
                        float red = std::min(d/200.0f, 1.0f);
                        drawCube((float)i + 0.5f, (float)j + 0.5f, (float)k + 0.5f,
                                0.85f, red, green, blue, alpha);
                    }
                }
            }
        }
    }

    // Draw grid velocity vectors (from actual grid data)
    if (state.showVelocityVectors) {
        for(int k = 1; k <= grid.N; k += state.vectorSkip) {
            for(int j = 1; j <= grid.N; j += state.vectorSkip) {
                for(int i = 1; i <= grid.N; i += state.vectorSkip) {
                    std::vector<float> vel = grid.getVelocityAtCellCenter(i, j, k);
                    float vx = vel[0];
                    float vy = vel[1];
                    float vz = vel[2];

                    float mag = std::sqrt(vx*vx + vy*vy + vz*vz);

                    // Only draw if velocity is significant
                    if (mag > state.minVelocityThreshold) {
                        float centerX = (float)i + 0.5f;
                        float centerY = (float)j + 0.5f;
                        float centerZ = (float)k + 0.5f;

                        // Scale velocity for visualization
                        float scale = state.vectorScale;
                        float endX = centerX + vx * scale;
                        float endY = centerY + vy * scale;
                        float endZ = centerZ + vz * scale;

                        // Color: cyan to yellow based on magnitude
                        float normalized = std::min(mag / 2.0f, 1.0f);
                        float r = normalized;
                        float g = 1.0f;
                        float b = 1.0f - normalized;

                        drawArrow(centerX, centerY, centerZ, endX, endY, endZ, r, g, b);
                    }
                }
            }
        }
    }

    // Draw interpolated velocity vectors (between grid points)
    if (state.showInterpolationVectors) {
        // First draw single interpolated vectors at specified offsets
        for(int k = 1; k <= grid.N; k += state.vectorSkip) {
            for(int j = 1; j <= grid.N; j += state.vectorSkip) {
                for(int i = 1; i <= grid.N; i += state.vectorSkip) {
                    // Test interpolation at offset position within the cell
                    float testX = (float)i + state.interpolationOffset;
                    float testY = (float)j + state.interpolationOffset;
                    float testZ = (float)k + state.interpolationOffset;

                    // Get velocity at cell center for comparison
                    std::vector<float> vel = grid.getVelocityAtCellCenter(i, j, k);
                    float centerMag = std::sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);

                    // Only show interpolation where there's significant velocity
                    if (centerMag > state.minVelocityThreshold * 0.5f) {
                        // Use the interpolation functions
                        float interpU = grid.interpolate_u(testX, testY, testZ);
                        float interpV = grid.interpolate_v(testX, testY, testZ);
                        float interpW = grid.interpolate_w(testX, testY, testZ);

                        float interpMag = std::sqrt(interpU*interpU + interpV*interpV + interpW*interpW);

                        if (interpMag > state.minVelocityThreshold * 0.5f) {
                            // Scale for visualization
                            float scale = state.vectorScale * 0.7f; // Slightly smaller
                            float endX = testX + interpU * scale;
                            float endY = testY + interpV * scale;
                            float endZ = testZ + interpW * scale;

                            // Color: magenta/purple to show it's interpolated
                            float normalized = std::min(interpMag / 2.0f, 1.0f);
                            float r = 1.0f;
                            float g = 0.3f;
                            float b = 1.0f - normalized * 0.5f;

                            drawArrow(testX, testY, testZ, endX, endY, endZ, r, g, b);
                        }
                    }
                }
            }
        }

        // Now draw subdivision interpolated vectors between grid points
        if (state.interpolationSubdivisions > 0) {
            // Calculate step size for subdivisions
            float step = 1.0f / (state.interpolationSubdivisions + 1);

            for(int k = 1; k < grid.N; k += state.vectorSkip) {
                for(int j = 1; j < grid.N; j += state.vectorSkip) {
                    for(int i = 1; i < grid.N; i += state.vectorSkip) {
                        // Create a mini grid between this cell and the next in each direction
                        for(int subZ = 0; subZ <= state.interpolationSubdivisions; subZ++) {
                            for(int subY = 0; subY <= state.interpolationSubdivisions; subY++) {
                                for(int subX = 0; subX <= state.interpolationSubdivisions; subX++) {
                                    // Position within the cell
                                    float subStepX = step * subX;
                                    float subStepY = step * subY;
                                    float subStepZ = step * subZ;

                                    // Skip if exactly at grid points (already shown above)
                                    if ((subX == 0 || subX == state.interpolationSubdivisions) &&
                                        (subY == 0 || subY == state.interpolationSubdivisions) &&
                                        (subZ == 0 || subZ == state.interpolationSubdivisions)) {
                                        continue;
                                    }

                                    float testX = (float)i + subStepX;
                                    float testY = (float)j + subStepY;
                                    float testZ = (float)k + subStepZ;

                                    // Interpolate velocity at this subdivision point
                                    float interpU = grid.interpolate_u(testX, testY, testZ);
                                    float interpV = grid.interpolate_v(testX, testY, testZ);
                                    float interpW = grid.interpolate_w(testX, testY, testZ);

                                    float interpMag = std::sqrt(interpU*interpU + interpV*interpV + interpW*interpW);

                                    if (interpMag > state.minVelocityThreshold * 0.3f) {
                                        // Scale for visualization
                                        float scale = state.vectorScale * 0.5f; // Smaller than regular vectors
                                        float endX = testX + interpU * scale;
                                        float endY = testY + interpV * scale;
                                        float endZ = testZ + interpW * scale;

                                        // Color: lighter magenta for subdivisions
                                        float normalized = std::min(interpMag / 2.0f, 1.0f);
                                        float r = 0.8f;
                                        float g = 0.6f;
                                        float b = 1.0f - normalized * 0.3f;

                                        drawArrow(testX, testY, testZ, endX, endY, endZ, r, g, b);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}



// --- Input Handling ---
void Input::handleCameraInput(CameraState& camera) {
    if (!ImGui::GetIO().WantCaptureMouse) {
        if (ImGui::IsMouseDown(ImGuiMouseButton_Right)) {
            ImVec2 mousePos = ImGui::GetMousePos();
            if (!camera.isDragging) {
                camera.isDragging = true;
                camera.lastMousePos = mousePos;
            }
            float dx = mousePos.x - camera.lastMousePos.x;
            float dy = mousePos.y - camera.lastMousePos.y;
            camera.rotY += dx * 0.5f;
            camera.rotX += dy * 0.5f;
            camera.lastMousePos = mousePos;
        } else {
            camera.isDragging = false;
        }

        camera.distance -= ImGui::GetIO().MouseWheel * 5.0f;
        if (camera.distance < 10.0f) camera.distance = 10.0f;
    }
}

// --- UI Implementation ---
void UI::setupImGui(GLFWwindow* window) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;

    if (!ImGui_ImplGlfw_InitForOpenGL(window, true)) {
        std::cerr << "Failed to initialize ImGui GLFW backend" << std::endl;
        return;
    }

    const char* glsl_version = "#version 130";
    if (!ImGui_ImplOpenGL3_Init(glsl_version)) {
        glsl_version = "#version 120";
        if (!ImGui_ImplOpenGL3_Init(glsl_version)) {
            std::cerr << "Failed to initialize ImGui OpenGL3 backend" << std::endl;
            return;
        }
    }

    ImGui::StyleColorsDark();
}

void UI::renderImGui(Grid3D& grid, SimulationState& state, CameraState& camera) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    ImGui::Begin("3D Fluid Controls");
    ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);

    // Simulation controls
    ImGui::Separator();
    ImGui::Text("Simulation:");
    ImGui::Checkbox("Running", &state.isRunning);
    ImGui::Checkbox("Show Divergence", &state.showDivergence);

    // Velocity visualization
    ImGui::Separator();
    ImGui::Text("Velocity Visualization:");
    ImGui::Checkbox("Show Grid Velocities", &state.showVelocityVectors);
    ImGui::Checkbox("Show Interpolated Velocities", &state.showInterpolationVectors);

    if (state.showVelocityVectors || state.showInterpolationVectors) {
        ImGui::SliderFloat("Vector Scale", &state.vectorScale, 0.1f, 10.0f);
        ImGui::SliderInt("Vector Skip", &state.vectorSkip, 1, 8);
        ImGui::SliderFloat("Min Velocity", &state.minVelocityThreshold, 0.0f, 1.0f);
    }

    if (state.showInterpolationVectors) {
        ImGui::SliderFloat("Interp Offset", &state.interpolationOffset, 0.1f, 0.9f);
        ImGui::SliderInt("Subdivisions", &state.interpolationSubdivisions, 0, 10);
        ImGui::Text("Cyan/Yellow = Grid data");
        ImGui::Text("Magenta = Interpolated (single)");
        ImGui::Text("Light Magenta = Interpolated (subdivided)");
    }

    // Velocity controls
    ImGui::Separator();
    ImGui::Text("Velocity Initialization (Staggered Grid):");

    if (ImGui::Button("Clear All Velocities")) {
        grid.clearAllVelocities();
    }

    if (ImGui::Button("Random Staggered")) {
        grid.initRandomStaggeredVelocities(0.5f);
    }

    // Density controls
    ImGui::Separator();
    ImGui::Text("Density Controls:");

    if (ImGui::Button("Add Random Density")) {
        int x = rand() % grid.N + 1;
        int y = rand() % grid.N + 1;
        int z = rand() % grid.N + 1;
        grid.add_density(x, y, z, 100.0f);
    }

    if (ImGui::Button("Add Random Velocity")) {
        int x = rand() % grid.N + 1;
        int y = rand() % grid.N + 1;
        int z = rand() % grid.N + 1;
        grid.add_velocity(x, y, z,
            (rand() / (float)RAND_MAX) * 2.0f - 1.0f,
            (rand() / (float)RAND_MAX) * 2.0f - 1.0f,
            (rand() / (float)RAND_MAX) * 2.0f - 1.0f);
    }

    if (ImGui::Button("Reset Everything")) {
        grid = Grid3D(Config::GRID_SIZE, state.diff, state.visc, state.dt);
    }

    // Simulation parameters
    ImGui::Separator();
    ImGui::Text("Parameters:");
    ImGui::SliderFloat("Viscosity", &state.visc, 0.0f, 0.001f, "%.5f");
    ImGui::SliderFloat("Diffusion", &state.diff, 0.0f, 0.001f, "%.5f");
    ImGui::SliderFloat("Time Step", &state.dt, 0.0f, 0.5f);

    // Camera controls
    ImGui::Separator();
    ImGui::Text("Camera:");
    ImGui::SliderFloat("Rotation X", &camera.rotX, -180.0f, 180.0f);
    ImGui::SliderFloat("Rotation Y", &camera.rotY, -180.0f, 180.0f);
    ImGui::SliderFloat("Distance", &camera.distance, 10.0f, 200.0f);

    // Divergence status
    ImGui::Separator();
    ImGui::Text("Status:");
    if (grid.checkDivergence(1e-2f)) {
        ImGui::TextColored(ImVec4(0, 1, 0, 1), "✓ Divergence: OK (incompressible)");
    } else {
        ImGui::TextColored(ImVec4(1, 0, 0, 1), "✗ Divergence: BAD (not incompressible)");
    }

    ImGui::Text("Grid Size: %d x %d x %d", grid.N, grid.N, grid.N);
    ImGui::Text("Using: Staggered MAC Grid");

    ImGui::End();

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

// --- App Implementation ---
GLFWwindow* App::initializeWindow() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return nullptr;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    GLFWwindow* window = glfwCreateWindow(
        Config::WINDOW_WIDTH,
        Config::WINDOW_HEIGHT,
        "3D Fluid Simulation (Staggered MAC Grid)",
        NULL, NULL
    );

    if (!window) {
        glfwTerminate();
        return nullptr;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        glfwDestroyWindow(window);
        glfwTerminate();
        return nullptr;
    }

    return window;
}

void App::initializeGraphics() {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
}

void App::mainLoop(GLFWwindow* window, Grid3D& grid, CameraState& camera, SimulationState& state) {
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // Handle input
        Input::handleCameraInput(camera);

        // Update simulation
        grid.dt = state.dt;
        grid.diff = state.diff;
        grid.visc = state.visc;

        if (state.isRunning) {
            grid.step();
        }

        // Get window dimensions
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);

        // Clear screen
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Setup view
        Graphics::setupPerspective(display_w, display_h);
        Graphics::setupCamera(camera);

        // Render 3D scene
        Graphics::renderScene(grid, state);

        // Render UI
        UI::renderImGui(grid, state, camera);

        // Swap buffers
        glfwSwapBuffers(window);
    }
}

void App::shutdown(GLFWwindow* window) {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
}

// --- Main Function ---
int main() {
    // Initialize window
    GLFWwindow* window = App::initializeWindow();
    if (!window) {
        return -1;
    }

    // Setup ImGui
    UI::setupImGui(window);

    // Setup graphics
    App::initializeGraphics();

    // Create simulation grid with staggered MAC grid
    Grid3D grid(Config::GRID_SIZE, 0.0001f, 0.0001f, 0.1f);

    // Initialize state
    CameraState camera;
    SimulationState state;

    std::cout << "Starting 3D Fluid Simulation with Staggered MAC Grid..." << std::endl;
    std::cout << "Grid configuration:" << std::endl;
    std::cout << "  - Pressure cells: " << Config::GRID_SIZE << "^3" << std::endl;
    std::cout << "  - u faces: " << (Config::GRID_SIZE+1) << " x " << (Config::GRID_SIZE+2) << " x " << (Config::GRID_SIZE+2) << std::endl;
    std::cout << "  - v faces: " << (Config::GRID_SIZE+2) << " x " << (Config::GRID_SIZE+1) << " x " << (Config::GRID_SIZE+2) << std::endl;
    std::cout << "  - w faces: " << (Config::GRID_SIZE+2) << " x " << (Config::GRID_SIZE+2) << " x " << (Config::GRID_SIZE+1) << std::endl;
    std::cout << "\nControls:" << std::endl;
    std::cout << "  - Right-click and drag to rotate camera" << std::endl;
    std::cout << "  - Mouse wheel to zoom" << std::endl;

    // Run main loop
    App::mainLoop(window, grid, camera, state);

    // Cleanup
    App::shutdown(window);

    std::cout << "Clean shutdown" << std::endl;
    return 0;
}