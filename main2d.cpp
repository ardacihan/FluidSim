#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "Grid2D.h"

// --- Configuration ---
int GRID_SIZE_2D = 64;
int WINDOW_WIDTH_2D = 1280;
int WINDOW_HEIGHT_2D = 720;

// Fixed color scaling values
float MAX_PRESSURE_2D = 1.0f;
float MAX_VELOCITY_2D = 2.0f;
float DENSITY_ALPHA_SCALE_2D = 0.8f;

// --- Global State ---
struct CameraState_2D {
    float zoom = 1.2f;
    float panX = 0.0f;
    float panY = 0.0f;
    ImVec2 lastMousePos;
    bool isPanning = false;
};

struct SimulationState_2D {
    float dt = 0.01f;
    float diff = 0.001f;
    float visc = 0.001f;
    bool isRunning = false;

    // Visualization modes
    enum VisualizationMode {
        DENSITY,
        PRESSURE,
        DYE
    };

    VisualizationMode visMode = DYE;

    enum FidelityLevel {
        LOW = 1,      // 1x1 (original)
        MEDIUM = 2,   // 2x2
        HIGH = 4,     // 4x4
        ULTRA = 8     // 8x8 (very high quality)
    };

    FidelityLevel fidelity = MEDIUM;

    // Velocity vector overlay
    bool showVelocityVectors = false;
    float vectorScale = 2.0f;
    int vectorSkip = 2;
    float minVelocityThreshold = 0.02f;

    // Mouse interaction
    ImVec2 lastInteractionPos;
    bool isAddingForce = false;
};

// --- Forward Declarations ---
namespace Graphics_2D {
    void drawSquare(float x, float y, float size, float r, float g, float b, float alpha);
    void drawWireframeBox(float width, float height);
    void drawArrow(float x1, float y1, float x2, float y2, float r, float g, float b);
    void setupOrthographic(int width, int height, const CameraState_2D& camera);
    void renderScene(Grid2D& grid, const SimulationState_2D& state);
    void drawCell(Grid2D& grid, const SimulationState_2D& state, int i, int j);
    void drawCellSupersampled(Grid2D& grid, const SimulationState_2D& state, int i, int j);
    void drawVelocityVectors(Grid2D& grid, const SimulationState_2D& state);
    std::vector<float> getDensityColor(float density);
    std::vector<float> getPressureColor(float pressure);
    std::vector<float> getDyeColor(float dyeValue);
}

namespace UI_2D {
    void setupImGUI(GLFWwindow* window);
    void renderImGUI(Grid2D& grid, SimulationState_2D& state, CameraState_2D& camera);
}

namespace Input {
    void handleCameraInput(CameraState_2D& camera);
    void handleMouseInteraction(GLFWwindow* window, Grid2D& grid, SimulationState_2D& state,
                                const CameraState_2D& camera, int width, int height);
}

namespace App_2D {
    GLFWwindow* initializeWindow();
    void initializeGraphics_2D();
    void mainLoop(GLFWwindow* window, Grid2D& grid, CameraState_2D& camera, SimulationState_2D& state);
    void shutdown(GLFWwindow* window);
}

// --- Graphics_2D Implementation ---
void Graphics_2D::drawSquare(float x, float y, float size, float r, float g, float b, float alpha) {
    glColor4f(r, g, b, alpha);

    float halfSize = size * 0.5f;
    glBegin(GL_QUADS);
    glVertex2f(x - halfSize, y - halfSize);
    glVertex2f(x + halfSize, y - halfSize);
    glVertex2f(x + halfSize, y + halfSize);
    glVertex2f(x - halfSize, y + halfSize);
    glEnd();
}

void Graphics_2D::drawWireframeBox(float width, float height) {
    glColor3f(1.0f, 1.0f, 1.0f);
    glLineWidth(2.0f);

    glBegin(GL_LINE_LOOP);
    glVertex2f(0, 0);
    glVertex2f(width, 0);
    glVertex2f(width, height);
    glVertex2f(0, height);
    glEnd();
}

void Graphics_2D::drawArrow(float x1, float y1, float x2, float y2, float r, float g, float b) {
    float dx = x2 - x1;
    float dy = y2 - y1;
    float length = std::sqrt(dx*dx + dy*dy);

    if (length < 0.001f) return;

    dx /= length;
    dy /= length;

    float arrowHeadLength = std::min(length * 0.3f, 0.5f);
    float arrowHeadWidth = arrowHeadLength * 0.5f;

    float baseX = x2 - dx * arrowHeadLength;
    float baseY = y2 - dy * arrowHeadLength;

    glColor3f(r, g, b);
    glLineWidth(2.0f);

    // Draw shaft
    glBegin(GL_LINES);
    glVertex2f(x1, y1);
    glVertex2f(baseX, baseY);
    glEnd();

    // Draw arrow head
    float perpX = -dy * arrowHeadWidth;
    float perpY = dx * arrowHeadWidth;

    glBegin(GL_TRIANGLES);
    glVertex2f(x2, y2);
    glVertex2f(baseX + perpX, baseY + perpY);
    glVertex2f(baseX - perpX, baseY - perpY);
    glEnd();
}

void Graphics_2D::setupOrthographic(int width, int height, const CameraState_2D& camera) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    float aspect = (float)width / (float)height;
    float viewWidth = GRID_SIZE_2D / camera.zoom;
    float viewHeight = viewWidth / aspect;

    float left = -viewWidth * 0.5f + camera.panX;
    float right = viewWidth * 0.5f + camera.panX;
    float bottom = viewHeight * 0.5f + camera.panY;
    float top = -viewHeight * 0.5f + camera.panY;

    glOrtho(left, right, bottom, top, -1.0f, 1.0f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}


std::vector<float> Graphics_2D::getDensityColor(float density) {
    std::vector<float> color(4, 0.0f); // RGBA

    if (density > 0.05f) {
        float alpha = std::min(density/100.0f, DENSITY_ALPHA_SCALE_2D);
        float normalized = std::min(density/200.0f, 1.0f);

        if (normalized < 0.5f) {
            float t = normalized * 2.0f;
            color = {0.0f, t, 1.0f, alpha};
        } else {
            float t = (normalized - 0.5f) * 2.0f;
            color = {t, 1.0f, 1.0f, alpha};
        }
    }

    return color;
}

std::vector<float> Graphics_2D::getPressureColor(float pressure) {
    std::vector<float> color(4, 0.0f); // RGBA

    float abs_p = std::abs(pressure);
    if (abs_p > 0.001f) {
        float normalized = std::min(abs_p / MAX_PRESSURE_2D, 1.0f);
        float alpha = 0.3f + 0.7f * normalized;
        color = {1.0f, 0.0f, 0.0f, alpha};
    }

    return color;
}

std::vector<float> Graphics_2D::getDyeColor(float dyeValue) {
    std::vector<float> color(4, 0.0f); // RGBA

    if (dyeValue > 0.5f) {
        float normalized = std::min(dyeValue / 100.0f, 1.0f);
        float alpha = std::min(normalized, 0.95f);

        // Beautiful color gradient: blue -> cyan -> green -> yellow -> red
        if (normalized < 0.25f) {
            float t = normalized * 4.0f;
            color = {0.0f, t, 1.0f, alpha};
        } else if (normalized < 0.5f) {
            float t = (normalized - 0.25f) * 4.0f;
            color = {0.0f, 1.0f, 1.0f - t, alpha};
        } else if (normalized < 0.75f) {
            float t = (normalized - 0.5f) * 4.0f;
            color = {t, 1.0f, 0.0f, alpha};
        } else {
            float t = (normalized - 0.75f) * 4.0f;
            color = {1.0f, 1.0f - t, 0.0f, alpha};
        }
    }

    return color;
}

void Graphics_2D::drawCell(Grid2D& grid, const SimulationState_2D& state, int i, int j) {
    std::vector<float> color;
    float value = 0.0f;

    switch (state.visMode) {
        case SimulationState_2D::DENSITY: {
            value = grid.dens[grid.P_IX(i, j)];
            color = getDensityColor(value);
            break;
        }

        case SimulationState_2D::PRESSURE: {
            value = grid.p[grid.P_IX(i, j)];
            color = getPressureColor(value);
            break;
        }

        case SimulationState_2D::DYE: {
            value = grid.dye[grid.P_IX(i, j)];
            color = getDyeColor(value);
            break;
        }
    }

    // Only draw if there's something to show
    if (color[3] > 0.01f) {
        drawSquare((float)i + 0.5f, (float)j + 0.5f, 0.9f,
                   color[0], color[1], color[2], color[3]);
    }
}

void Graphics_2D::drawCellSupersampled(Grid2D& grid, const SimulationState_2D& state, int i, int j) {
    int subdiv = (int)state.fidelity;

    // Early out for low fidelity (original behavior)
    if (subdiv == 1) {
        drawCell(grid, state, i, j);
        return;
    }

    // High fidelity: supersampling
    const float cellSize = 1.0f;
    const float subSize = cellSize / subdiv;
    const float drawSize = subSize * 0.9f;

    float baseX = (float)i + 0.5f;
    float baseY = (float)j + 0.5f;

    // Use immediate mode for simplicity
    for (int subY = 0; subY < subdiv; subY++) {
        float subCenterY = baseY - 0.5f + (subY + 0.5f) * subSize;

        for (int subX = 0; subX < subdiv; subX++) {
            float subCenterX = baseX - 0.5f + (subX + 0.5f) * subSize;

            // Get interpolated value at this sub-pixel position
            float interpolatedValue = 0.0f;
            std::vector<float> color(4, 0.0f);

            switch (state.visMode) {
                case SimulationState_2D::DENSITY: {
                    interpolatedValue = grid.interpolate_density(subCenterX, subCenterY);
                    color = getDensityColor(interpolatedValue);
                    break;
                }

                case SimulationState_2D::PRESSURE: {
                    // Need to add pressure interpolation to Grid2D or use cell value
                    // For now, use cell value
                    interpolatedValue = grid.p[grid.P_IX(i, j)];
                    color = getPressureColor(interpolatedValue);
                    break;
                }

                case SimulationState_2D::DYE: {
                    interpolatedValue = grid.interpolate_dye(subCenterX, subCenterY);
                    color = getDyeColor(interpolatedValue);
                    break;
                }
            }

            // Only draw if there's something to show
            if (color[3] > 0.01f) {
                drawSquare(subCenterX, subCenterY, drawSize,
                           color[0], color[1], color[2], color[3]);
            }
        }
    }
}

void Graphics_2D::drawVelocityVectors(Grid2D& grid, const SimulationState_2D& state) {
    for(int j = 1; j <= grid.N; j += state.vectorSkip) {
        for(int i = 1; i <= grid.N; i += state.vectorSkip) {
            std::vector<float> vel = grid.getVelocityAtCellCenter(i, j);
            float vx = vel[0];
            float vy = vel[1];
            float mag = std::sqrt(vx*vx + vy*vy);

            if (mag > state.minVelocityThreshold) {
                float centerX = (float)i + 0.5f;
                float centerY = (float)j + 0.5f;
                float scale = state.vectorScale;
                float endX = centerX + vx * scale;
                float endY = centerY + vy * scale;

                float normalized = std::min(mag / MAX_VELOCITY_2D, 1.0f);
                float r = normalized;
                float g = 0.7f;
                float b = 1.0f - normalized;

                drawArrow(centerX, centerY, endX, endY, r, g, b);
            }
        }
    }
}

void Graphics_2D::renderScene(Grid2D& grid, const SimulationState_2D& state) {
    // Draw grid boundary
    drawWireframeBox((float)GRID_SIZE_2D, (float)GRID_SIZE_2D);

    // Draw cells with supersampling if enabled
    for(int j = 1; j <= grid.N; j++) {
        for(int i = 1; i <= grid.N; i++) {
            if (state.fidelity == SimulationState_2D::LOW) {
                drawCell(grid, state, i, j);
            } else {
                drawCellSupersampled(grid, state, i, j);
            }
        }
    }

    // Draw velocity vectors if enabled
    if (state.showVelocityVectors) {
        drawVelocityVectors(grid, state);
    }
}

// --- Input Handling ---
void Input::handleCameraInput(CameraState_2D& camera) {
    if (!ImGui::GetIO().WantCaptureMouse) {
        // Pan with middle mouse button
        if (ImGui::IsMouseDown(ImGuiMouseButton_Middle)) {
            ImVec2 mousePos = ImGui::GetMousePos();
            if (!camera.isPanning) {
                camera.isPanning = true;
                camera.lastMousePos = mousePos;
            }
            float dx = mousePos.x - camera.lastMousePos.x;
            float dy = mousePos.y - camera.lastMousePos.y;

            camera.panX -= dx * (GRID_SIZE_2D / camera.zoom) * 0.01f;
            camera.panY += dy * (GRID_SIZE_2D / camera.zoom) * 0.01f;
            camera.lastMousePos = mousePos;
        } else {
            camera.isPanning = false;
        }

        float wheel = ImGui::GetIO().MouseWheel;
        if (wheel != 0) {
            camera.zoom *= (1.0f + wheel * 0.1f);
            camera.zoom = std::max(0.5f, std::min(camera.zoom, 20.0f));
        }
    }
}

void Input::handleMouseInteraction(GLFWwindow* window, Grid2D& grid, SimulationState_2D& state,
                                   const CameraState_2D& camera, int width, int height) {
    if (ImGui::GetIO().WantCaptureMouse) return;

    // Check if CTRL + Left Mouse Button is pressed
    bool ctrlPressed = glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
                       glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;

    if (ctrlPressed && ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
        ImVec2 mousePos = ImGui::GetMousePos();

        // Convert screen coordinates to grid coordinates
        float aspect = (float)width / (float)height;
        float viewWidth = GRID_SIZE_2D / camera.zoom;
        float viewHeight = viewWidth / aspect;

        // Calculate normalized coordinates (0 to 1)
        float normX = mousePos.x / width;
        float normY = mousePos.y / height;

        // Map to grid coordinates
        float left = -viewWidth * 0.5f + camera.panX;
        float top = -viewHeight * 0.5f + camera.panY;

        float gridPosX = left + normX * viewWidth;
        float gridPosY = top + normY * viewHeight;

        int gridX = (int)gridPosX;
        int gridY = (int)gridPosY;

        // Boundary check
        if (gridX >= 1 && gridX <= grid.N && gridY >= 1 && gridY <= grid.N) {
            if (!state.isAddingForce) {
                state.isAddingForce = true;
                state.lastInteractionPos = mousePos;
            }

            // Calculate mouse velocity
            float dx = mousePos.x - state.lastInteractionPos.x;
            float dy = mousePos.y - state.lastInteractionPos.y;

            // Convert screen velocity to world velocity
            float velX = dx * (GRID_SIZE_2D) * 0.01f;
            float velY = dy * (GRID_SIZE_2D) * 0.01f;

            // Add velocity and density
            grid.add_velocity(gridX, gridY, velX, velY);
            grid.add_density(gridX, gridY, 300.0f);
            grid.add_dye(gridX, gridY, 300.0f);

            state.lastInteractionPos = mousePos;
        }
    } else {
        state.isAddingForce = false;
    }
}

// --- UI_2D Implementation ---
void UI_2D::setupImGUI(GLFWwindow* window) {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    (void)io;

    ImGui_ImplGlfw_InitForOpenGL(window, true);
    const char* glsl_version = "#version 130";
    ImGui_ImplOpenGL3_Init(glsl_version);

    ImGui::StyleColorsDark();
}

void UI_2D::renderImGUI(Grid2D& grid, SimulationState_2D& state, CameraState_2D& camera) {
    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();

    ImGui::Begin("2D Fluid Controls");

    // FPS and stats
    ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
    ImGui::Text("Grid: %dx%d cells", grid.N, grid.N);

    // Simulation controls
    ImGui::Separator();
    ImGui::Text("Simulation Control:");
    ImGui::Checkbox("Running", &state.isRunning);
    ImGui::SameLine();
    if (ImGui::Button("Single Step")) {
        grid.step();
    }

    // Timestep parameters
    ImGui::Separator();
    ImGui::Text("Simulation Parameters:");
    ImGui::SliderFloat("Time Step (dt)", &state.dt, 0.001f, 0.1f, "%.4f");
    ImGui::SliderFloat("Viscosity", &state.visc, 0.0f, 0.1f, "%.4f");
    ImGui::SliderFloat("Density Diffusion", &state.diff, 0.0f, 0.1f, "%.4f");

    ImGui::Separator();
    ImGui::Text("Visualization Quality:");
    if (ImGui::RadioButton("Low (1x1)", state.fidelity == SimulationState_2D::LOW)) {
        state.fidelity = SimulationState_2D::LOW;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("Medium (2x2)", state.fidelity == SimulationState_2D::MEDIUM)) {
        state.fidelity = SimulationState_2D::MEDIUM;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("High (4x4)", state.fidelity == SimulationState_2D::HIGH)) {
        state.fidelity = SimulationState_2D::HIGH;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("Ultra (8x8)", state.fidelity == SimulationState_2D::ULTRA)) {
        state.fidelity = SimulationState_2D::ULTRA;
    }

    // Visualization mode selection
    ImGui::Separator();
    ImGui::Text("Visualization:");

    if (ImGui::RadioButton("Density Field", state.visMode == SimulationState_2D::DENSITY)) {
        state.visMode = SimulationState_2D::DENSITY;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("Pressure Field", state.visMode == SimulationState_2D::PRESSURE)) {
        state.visMode = SimulationState_2D::PRESSURE;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("Dye", state.visMode == SimulationState_2D::DYE)) {
        state.visMode = SimulationState_2D::DYE;
    }

    // Show color scheme info
    ImGui::Separator();
    ImGui::Text("Color Scheme:");
    if (state.visMode == SimulationState_2D::DENSITY) {
        ImGui::TextColored(ImVec4(0, 0, 1, 1), "Blue: Low density");
        ImGui::TextColored(ImVec4(0, 1, 1, 1), "Cyan: Medium density");
        ImGui::TextColored(ImVec4(1, 1, 1, 1), "White: High density");
    } else if (state.visMode == SimulationState_2D::PRESSURE) {
        ImGui::TextColored(ImVec4(1, 0, 0, 0.3f), "Transparent Red: Low pressure");
        ImGui::TextColored(ImVec4(1, 0, 0, 1.0f), "Opaque Red: High pressure");
    } else {
        ImGui::TextColored(ImVec4(0, 0, 1, 1), "Blue -> Cyan -> Green -> Yellow -> Red");
        ImGui::Text("(Based on dye concentration)");
    }

    // Vector visualization
    ImGui::Separator();
    ImGui::Text("Overlay:");
    ImGui::Checkbox("Show Velocity Vectors", &state.showVelocityVectors);

    if (state.showVelocityVectors) {
        ImGui::SliderFloat("Vector Scale", &state.vectorScale, 0.1f, 100.0f);
        ImGui::SliderInt("Vector Skip", &state.vectorSkip, 1, 8);
        ImGui::SliderFloat("Min Velocity", &state.minVelocityThreshold, 0.0f, 1.0f);
    }

    // Mouse interaction
    ImGui::Separator();
    ImGui::Text("Mouse Interaction:");
    ImGui::TextColored(ImVec4(0.5f, 1.0f, 0.5f, 1), "CTRL + Left Click: Add velocity & density");
    ImGui::TextColored(ImVec4(0.5f, 1.0f, 0.5f, 1), "Middle Mouse: Pan view");
    ImGui::TextColored(ImVec4(0.5f, 1.0f, 0.5f, 1), "Mouse Wheel: Zoom");

    // Initialization controls
    ImGui::Separator();
    ImGui::Text("Initialization:");

    if (ImGui::Button("Clear All")) {
        grid.clearAllVelocities();
        std::fill(grid.dens.begin(), grid.dens.end(), 0.0f);
        std::fill(grid.p.begin(), grid.p.end(), 0.0f);
        grid.clearAllDye();
    }
    ImGui::SameLine();
    if (ImGui::Button("Random Velocities")) {
        grid.initRandomStaggeredVelocities(1.0f);
    }

    // Camera controls
    ImGui::Separator();
    ImGui::Text("Camera:");
    ImGui::SliderFloat("Zoom", &camera.zoom, 0.5f, 20.0f);
    ImGui::SliderFloat("Pan X", &camera.panX, -GRID_SIZE_2D*2.0f, GRID_SIZE_2D*2.0f);
    ImGui::SliderFloat("Pan Y", &camera.panY, -GRID_SIZE_2D*2.0f, GRID_SIZE_2D*2.0f);

    if (ImGui::Button("Reset Camera")) {
        camera.zoom = 1.2f;
        camera.panX = GRID_SIZE_2D * 0.5f;
        camera.panY = GRID_SIZE_2D * 0.5f;
    }

    // Statistics
    ImGui::Separator();
    ImGui::Text("Statistics:");

    float maxVel = 0.0f;
    float maxDens = 0.0f;
    float maxPress = 0.0f;
    float totalDye = 0.0f;

    for (int i = 0; i < grid.u.size(); i++) {
        if (fabs(grid.u[i]) > maxVel) maxVel = fabs(grid.u[i]);
    }
    for (int i = 0; i < grid.v.size(); i++) {
        if (fabs(grid.v[i]) > maxVel) maxVel = fabs(grid.v[i]);
    }

    for (int i = 0; i < grid.dens.size(); i++) {
        if (grid.dens[i] > maxDens) maxDens = grid.dens[i];
    }

    for (int i = 0; i < grid.p.size(); i++) {
        if (fabs(grid.p[i]) > maxPress) maxPress = fabs(grid.p[i]);
    }

    for (int i = 0; i < grid.dye.size(); i++) {
        totalDye += grid.dye[i];
    }

    ImGui::Text("Max Velocity: %.4f", maxVel);
    ImGui::Text("Max Density: %.2f", maxDens);
    ImGui::Text("Max Pressure: %.4f", maxPress);
    ImGui::Text("Total Dye: %.2f", totalDye);

    // Reset button
    ImGui::Separator();
    if (ImGui::Button("Reset Everything")) {
        grid = Grid2D(GRID_SIZE_2D, state.diff, state.visc, state.dt);
        state.isRunning = false;
        state.visMode = SimulationState_2D::DYE;
    }

    ImGui::End();

    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
}

// --- App Implementation ---
GLFWwindow* App_2D::initializeWindow() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return nullptr;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);

    GLFWwindow* window = glfwCreateWindow(
        WINDOW_WIDTH_2D,
        WINDOW_HEIGHT_2D,
        "2D Fluid Simulation with Supersampling",
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
        std::cerr << "Failed to initialize GLEW" << std::endl;
        glfwDestroyWindow(window);
        glfwTerminate();
        return nullptr;
    }

    return window;
}

void App_2D::initializeGraphics_2D() {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
}

void App_2D::mainLoop(GLFWwindow* window, Grid2D& grid, CameraState_2D& camera, SimulationState_2D& state) {
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // Get window dimensions
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);

        // Handle input
        Input::handleCameraInput(camera);
        Input::handleMouseInteraction(window, grid, state, camera, display_w, display_h);

        // Update simulation
        grid.dt = state.dt;
        grid.diff = state.diff;
        grid.visc = state.visc;

        if (state.isRunning) {
            grid.step();
        }

        // Clear screen
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // Setup view BEFORE rendering
        Graphics_2D::setupOrthographic(display_w, display_h, camera);

        // Render 2D scene
        Graphics_2D::renderScene(grid, state);

        // Render UI (this should be last)
        UI_2D::renderImGUI(grid, state, camera);

        // Swap buffers
        glfwSwapBuffers(window);
    }
}

void App_2D::shutdown(GLFWwindow* window) {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
}

// --- Main Function ---
int main() {
    // Initialize window
    GLFWwindow* window = App_2D::initializeWindow();
    if (!window) {
        return -1;
    }

    // Setup ImGUI
    UI_2D::setupImGUI(window);

    // Setup Graphics
    App_2D::initializeGraphics_2D();

    // Create simulation grid
    Grid2D grid(GRID_SIZE_2D, 0.001f, 0.001f, 0.01f);

    // Initialize state
    CameraState_2D camera;
    camera.panX = GRID_SIZE_2D * 0.5f;
    camera.panY = GRID_SIZE_2D * 0.5f;

    SimulationState_2D state;

    // Run main loop
    App_2D::mainLoop(window, grid, camera, state);

    // Cleanup
    App_2D::shutdown(window);
    return 0;
}