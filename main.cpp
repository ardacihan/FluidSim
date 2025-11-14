#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Grid.h"

// =============================================================================
// CONFIGURATION - Adjust these values to customize the simulation
// =============================================================================

// Grid dimensions
const int GRID_COLS = 192;
const int GRID_ROWS = 108;

// Window settings - High resolution fullscreen
const int WINDOW_WIDTH = 1920;
const int WINDOW_HEIGHT = 1080;
const char* WINDOW_TITLE = "Fluid Simulation - Staggered Grid";

// Grid colors (R, G, B, A) - values from 0-255
const ImU32 BACKGROUND_TOP = IM_COL32(245, 247, 250, 255);
const ImU32 BACKGROUND_BOTTOM = IM_COL32(235, 237, 242, 255);
const ImU32 INNER_LINE_COLOR = IM_COL32(180, 185, 195, 255);
const ImU32 WALL_COLOR = IM_COL32(35, 35, 40, 255);
const float INNER_LINE_THICKNESS = 1.5f;
const float WALL_THICKNESS = 4.0f;

// Arrow visualization
const float ARROW_MAX_LENGTH_RATIO = 0.8f;  // Max arrow length as fraction of cell size
const float ARROW_BASE_HEAD_SIZE = 8.0f;    // Base size of arrow heads
const float ARROW_MIN_WIDTH = 1.5f;         // Minimum arrow line width
const float ARROW_MAX_WIDTH = 3.5f;         // Maximum arrow line width
const float FLOW_THRESHOLD = 0.01f;         // Minimum flow to display arrow

// Arrow colors (low to high intensity)
const int ARROW_COLOR_LOW_R = 50;
const int ARROW_COLOR_LOW_G = 120;
const int ARROW_COLOR_LOW_B = 220;
const int ARROW_COLOR_HIGH_R = 220;
const int ARROW_COLOR_HIGH_G = 80;
const int ARROW_COLOR_HIGH_B = 80;

// UI settings
const float UI_FONT_SCALE = 1.3f;
const float UI_WINDOW_ROUNDING = 8.0f;
const float UI_FRAME_ROUNDING = 4.0f;

// Mouse interaction
const float MOUSE_DENSITY_AMOUNT = 10.0f;     // How much density to add per click
const float MOUSE_VELOCITY_STRENGTH = 1.0f; // How much velocity to add per drag
const float MOUSE_RADIUS = 1.0f;             // Radius of effect in cells

// Simulation parameters
const float DIFFUSION_RATE = 0.0001f;
const float VISCOSITY = 0.0001f;
const float TIME_STEP = 0.05f;

// =============================================================================
// Helper Functions
// =============================================================================

void initializeImGuiStyle() {
    ImGuiStyle& style = ImGui::GetStyle();
    style.WindowRounding = UI_WINDOW_ROUNDING;
    style.FrameRounding = UI_FRAME_ROUNDING;
    style.GrabRounding = UI_FRAME_ROUNDING;
    style.Colors[ImGuiCol_WindowBg] = ImVec4(0.12f, 0.12f, 0.15f, 0.95f);
    style.Colors[ImGuiCol_Button] = ImVec4(0.2f, 0.4f, 0.7f, 1.0f);
    style.Colors[ImGuiCol_ButtonHovered] = ImVec4(0.3f, 0.5f, 0.8f, 1.0f);
    style.Colors[ImGuiCol_ButtonActive] = ImVec4(0.15f, 0.35f, 0.65f, 1.0f);
}

// Calculate optimal cell size to fit grid perfectly in window
float calculateCellSize(int windowWidth, int windowHeight, int gridCols, int gridRows) {
    // Grid includes 2 extra cells for walls on each side
    int totalCols = gridCols + 2;
    int totalRows = gridRows + 2;

    // Calculate maximum cell size that fits both dimensions
    float cellSizeByWidth = (float)windowWidth / totalCols;
    float cellSizeByHeight = (float)windowHeight / totalRows;

    // Use the smaller of the two to ensure it fits
    return std::min(cellSizeByWidth, cellSizeByHeight);
}

void drawCellBackgrounds(ImDrawList* draw_list, const Grid& grid, float cellSize,
                        float startX, float startY) {
    // Find max density for normalization
    float maxDensity = 0.0f;
    for (int i = 1; i <= grid.CellCountX; ++i) {
        for (int j = 1; j <= grid.CellCountY; ++j) {
            maxDensity = std::max(maxDensity, grid.d[i][j]);
        }
    }

    // Avoid division by zero
    if (maxDensity < 0.001f) maxDensity = 1.0f;

    // Draw density visualization: black (low) to white (high)
    for (int i = 1; i <= grid.CellCountX; ++i) {
        for (int j = 1; j <= grid.CellCountY; ++j) {
            float x1 = startX + i * cellSize;
            float y1 = startY + j * cellSize;
            float x2 = x1 + cellSize;
            float y2 = y1 + cellSize;

            // Normalize density to [0, 1]
            float normalizedDensity = grid.d[i][j] / maxDensity;
            normalizedDensity = std::clamp(normalizedDensity, 0.0f, 1.0f);

            // Map to grayscale: 0 (black) to 255 (white)
            int gray = (int)(normalizedDensity * 255);
            ImU32 color = IM_COL32(gray, gray, gray, 255);

            draw_list->AddRectFilled(ImVec2(x1, y1), ImVec2(x2, y2), color);
        }
    }
}

void drawGridLines(ImDrawList* draw_list, const Grid& grid, float cellSize,
                   float startX, float startY, float gridWidth, float gridHeight) {
    // Vertical lines
    for (int i = 0; i <= grid.CellCountX + 2; ++i) {
        float x = startX + i * cellSize;
        bool isWall = (i == 0 || i == grid.CellCountX + 2);
        float thickness = isWall ? WALL_THICKNESS : INNER_LINE_THICKNESS;
        ImU32 color = isWall ? WALL_COLOR : INNER_LINE_COLOR;

        draw_list->AddLine(
            ImVec2(x, startY),
            ImVec2(x, startY + gridHeight),
            color, thickness
        );
    }

    // Horizontal lines
    for (int i = 0; i <= grid.CellCountY + 2; ++i) {
        float y = startY + i * cellSize;
        bool isWall = (i == 0 || i == grid.CellCountY + 2);
        float thickness = isWall ? WALL_THICKNESS : INNER_LINE_THICKNESS;
        ImU32 color = isWall ? WALL_COLOR : INNER_LINE_COLOR;

        draw_list->AddLine(
            ImVec2(startX, y),
            ImVec2(startX + gridWidth, y),
            color, thickness
        );
    }
}

void drawFlowVectors(ImDrawList* draw_list, const Grid& grid, float cellSize,
                     float startX, float startY) {
    // Draw arrows at cell centers using interpolated velocities
    for (int i = 1; i <= grid.CellCountX; ++i) {
        for (int j = 1; j <= grid.CellCountY; ++j) {
            // Calculate cell center position
            float cx = startX + (i + 0.5f) * cellSize;
            float cy = startY + (j + 0.5f) * cellSize;

            // Interpolate velocities to cell center
            float u_center = 0.5f * (grid.u[i-1][j] + grid.u[i][j]);
            float v_center = 0.5f * (grid.v[i][j-1] + grid.v[i][j]);

            // Calculate magnitude
            float magnitude = std::sqrt(u_center * u_center + v_center * v_center);

            if (magnitude > FLOW_THRESHOLD) {
                // Normalize direction
                float dirX = u_center / magnitude;
                float dirY = v_center / magnitude;

                // Scale arrow length based on magnitude
                float maxLength = ARROW_MAX_LENGTH_RATIO * cellSize;
                float arrowLength = std::min(magnitude, 1.0f) * maxLength;

                // Calculate arrow end point
                float ex = cx + dirX * arrowLength;
                float ey = cy + dirY * arrowLength;

                // Arrow head size scales with magnitude
                float arrowHeadSize = ARROW_BASE_HEAD_SIZE * (0.5f + 0.5f * std::min(magnitude, 1.0f));

                // Color gradient from low (blue) to high (red) intensity
                float intensity = std::min(magnitude, 1.0f);
                int r = (int)(ARROW_COLOR_LOW_R + intensity * (ARROW_COLOR_HIGH_R - ARROW_COLOR_LOW_R));
                int g = (int)(ARROW_COLOR_LOW_G + intensity * (ARROW_COLOR_HIGH_G - ARROW_COLOR_LOW_G));
                int b = (int)(ARROW_COLOR_LOW_B + intensity * (ARROW_COLOR_HIGH_B - ARROW_COLOR_LOW_B));
                ImU32 color = IM_COL32(r, g, b, (int)(intensity * 200 + 55));

                // Line width scales with magnitude
                float lineWidth = ARROW_MIN_WIDTH + intensity * (ARROW_MAX_WIDTH - ARROW_MIN_WIDTH);

                // Draw arrow shaft
                draw_list->AddLine(ImVec2(cx, cy), ImVec2(ex, ey), color, lineWidth);

                // Calculate arrow head points
                float perpX = -dirY;
                float perpY = dirX;

                ImVec2 tip(ex, ey);
                ImVec2 wing1(ex - dirX * arrowHeadSize + perpX * arrowHeadSize * 0.5f,
                            ey - dirY * arrowHeadSize + perpY * arrowHeadSize * 0.5f);
                ImVec2 wing2(ex - dirX * arrowHeadSize - perpX * arrowHeadSize * 0.5f,
                            ey - dirY * arrowHeadSize - perpY * arrowHeadSize * 0.5f);

                draw_list->AddTriangleFilled(tip, wing1, wing2, color);
            }
        }
    }
}

void addDensityAtMouse(Grid& grid, float mouseX, float mouseY, float startX, float startY, float cellSize) {
    // Convert mouse position to grid coordinates
    float gridX = (mouseX - startX) / cellSize;
    float gridY = (mouseY - startY) / cellSize;

    // Check if mouse is within grid bounds (excluding walls)
    if (gridX >= 1.0f && gridX <= grid.CellCountX + 1.0f &&
        gridY >= 1.0f && gridY <= grid.CellCountY + 1.0f) {

        int centerI = static_cast<int>(gridX);
        int centerJ = static_cast<int>(gridY);

        // Add density in a circular area around the mouse
        for (int i = std::max(1, centerI - static_cast<int>(MOUSE_RADIUS));
             i <= std::min(grid.CellCountX, centerI + static_cast<int>(MOUSE_RADIUS));
             ++i) {
            for (int j = std::max(1, centerJ - static_cast<int>(MOUSE_RADIUS));
                 j <= std::min(grid.CellCountY, centerJ + static_cast<int>(MOUSE_RADIUS));
                 ++j) {

                // Calculate distance from mouse center
                float dist = std::sqrt((i - gridX) * (i - gridX) + (j - gridY) * (j - gridY));

                if (dist <= MOUSE_RADIUS) {
                    // Add more density closer to center (inverse square falloff)
                    float falloff = 1.0f - (dist / MOUSE_RADIUS);
                    grid.d[i][j] += MOUSE_DENSITY_AMOUNT * falloff * falloff;
                }
            }
        }
    }
}

void addVelocityAtMouse(Grid& grid, float mouseX, float mouseY, float velX, float velY,
                       float startX, float startY, float cellSize) {
    // Convert mouse position to grid coordinates
    float gridX = (mouseX - startX) / cellSize;
    float gridY = (mouseY - startY) / cellSize;

    // Check if mouse is within grid bounds (excluding walls)
    if (gridX >= 1.0f && gridX <= grid.CellCountX + 1.0f &&
        gridY >= 1.0f && gridY <= grid.CellCountY + 1.0f) {

        int centerI = static_cast<int>(gridX);
        int centerJ = static_cast<int>(gridY);

        // Add velocity in a circular area around the mouse
        for (int i = std::max(1, centerI - static_cast<int>(MOUSE_RADIUS));
             i <= std::min(grid.CellCountX, centerI + static_cast<int>(MOUSE_RADIUS));
             ++i) {
            for (int j = std::max(1, centerJ - static_cast<int>(MOUSE_RADIUS));
                 j <= std::min(grid.CellCountY, centerJ + static_cast<int>(MOUSE_RADIUS));
                 ++j) {

                // Calculate distance from mouse center
                float dist = std::sqrt((i - gridX) * (i - gridX) + (j - gridY) * (j - gridY));

                if (dist <= MOUSE_RADIUS) {
                    // Add more velocity closer to center (inverse square falloff)
                    float falloff = 1.0f - (dist / MOUSE_RADIUS);
                    float strength = falloff * falloff;

                    // Apply velocity to both u and v components
                    grid.u[i][j] += velX * strength;
                    grid.v[i][j] += velY * strength;
                }
            }
        }
    }
}

void renderControlPanel(Grid& grid, bool& simulationRunning, float& timeStep, float& diffusionRate, float& viscosity,
                       bool& showArrows, const ImGuiIO& io, float startX, float startY, float cellSize) {
    ImGui::SetNextWindowPos(ImVec2(30, 30), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(400, 450), ImGuiCond_FirstUseEver);
    ImGui::Begin("Fluid Simulation Controls", nullptr, ImGuiWindowFlags_NoCollapse);

    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.9f, 0.9f, 0.95f, 1.0f));

    ImGui::Text("Grid Configuration: %dx%d cells", grid.CellCountX, grid.CellCountY);
    ImGui::Text("Cell Size: %.1f px (auto-fitted)", grid.cellDisplaySize);
    ImGui::Separator();
    ImGui::Spacing();

    ImGui::Text("Visualization:");
    ImGui::Checkbox("Show Velocity Arrows", &showArrows);
    ImGui::Indent();
    ImGui::Text("Density: Black (low) → White (high)");
    ImGui::Text("Arrows: Blue (slow) → Red (fast)");
    ImGui::Unindent();

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    ImGui::Text("Simulation Parameters:");
    ImGui::SliderFloat("Time Step", &timeStep, 0.01f, 0.5f, "%.3f");
    ImGui::SliderFloat("Diffusion Rate", &diffusionRate, 0.0f, 0.1f, "%.4f");
    ImGui::SliderFloat("Viscosity", &viscosity, 0.0f, 0.1f, "%.4f");

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    ImGui::Text("Controls:");

    // Simulation control buttons
    if (ImGui::Button(simulationRunning ? "Pause Simulation" : "Start Simulation", ImVec2(180, 35))) {
        simulationRunning = !simulationRunning;
    }
    ImGui::SameLine();
    if (ImGui::Button("Step Once", ImVec2(180, 35))) {
        // Single simulation step
        grid.vel_step(viscosity, timeStep);
        grid.density_step(diffusionRate, timeStep);
    }

    if (ImGui::Button("Randomize Density", ImVec2(180, 35))) {
        grid.initializeDensity();
    }
    ImGui::SameLine();
    if (ImGui::Button("Clear Density", ImVec2(180, 35))) {
        grid.clearDensity();
    }

    ImGui::Spacing();
    if (ImGui::Button("Randomize Velocities", ImVec2(180, 35))) {
        grid.initializeRandomVelocities();
    }
    ImGui::SameLine();
    if (ImGui::Button("Clear Velocities", ImVec2(180, 35))) {
        grid.resetVelocities();
    }

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    ImGui::Text("Mouse Interaction:");
    ImGui::Text("Left Click: Add density");
    ImGui::Text("Right Click + Drag: Add velocity");
    ImGui::Text("Grid Area: (%.0f, %.0f) to (%.0f, %.0f)",
                startX, startY, startX + (grid.CellCountX + 2) * cellSize, startY + (grid.CellCountY + 2) * cellSize);

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    ImGui::Text("Performance");
    ImGui::Text("FPS: %.1f", io.Framerate);

    ImGui::PopStyleColor();
    ImGui::End();
}

// =============================================================================
// Main Program
// =============================================================================

int main() {
    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, WINDOW_TITLE, nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    io.FontGlobalScale = UI_FONT_SCALE;

    ImGui::StyleColorsDark();
    initializeImGuiStyle();

    if (!ImGui_ImplGlfw_InitForOpenGL(window, true)) {
        std::cerr << "Failed to initialize ImGui GLFW backend" << std::endl;
        return -1;
    }

    if (!ImGui_ImplOpenGL3_Init("#version 330")) {
        std::cerr << "Failed to initialize ImGui OpenGL backend" << std::endl;
        return -1;
    }

    // Initialize simulation grid
    Grid grid(GRID_COLS, GRID_ROWS, 1.0f);

    // Calculate cell size to fit window perfectly
    float cellSize = calculateCellSize(WINDOW_WIDTH, WINDOW_HEIGHT, GRID_COLS, GRID_ROWS);
    grid.cellDisplaySize = cellSize;

    // Simulation state
    bool simulationRunning = false;
    bool showArrows = true;

    // Mouse drag tracking for velocity
    bool isRightDragging = false;
    ImVec2 lastMousePos = ImVec2(0, 0);

    // Simulation parameters with default values
    float timeStep = TIME_STEP;
    float diffusionRate = DIFFUSION_RATE;
    float viscosity = VISCOSITY;

    // Main loop
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Get window size and setup viewport
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.95f, 0.96f, 0.97f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        ImDrawList* draw_list = ImGui::GetBackgroundDrawList();

        // Calculate grid dimensions and centering
        float gridWidth = (GRID_COLS + 2) * cellSize;
        float gridHeight = (GRID_ROWS + 2) * cellSize;
        float startX = (display_w - gridWidth) / 2.0f;
        float startY = (display_h - gridHeight) / 2.0f;

        // Handle mouse input
        ImVec2 currentMousePos = ImGui::GetMousePos();

        // Left click - add density
        if (ImGui::IsMouseDown(ImGuiMouseButton_Left) && !ImGui::GetIO().WantCaptureMouse) {
            addDensityAtMouse(grid, currentMousePos.x, currentMousePos.y, startX, startY, cellSize);
        }

        // Right click drag - add velocity
        if (ImGui::IsMouseDown(ImGuiMouseButton_Right) && !ImGui::GetIO().WantCaptureMouse) {
            if (!isRightDragging) {
                // Start dragging
                isRightDragging = true;
                lastMousePos = currentMousePos;
            } else {
                // Calculate drag velocity (direction and magnitude)
                float dragX = currentMousePos.x - lastMousePos.x;
                float dragY = currentMousePos.y - lastMousePos.y;

                // Normalize and scale by strength
                float length = std::sqrt(dragX * dragX + dragY * dragY);
                if (length > 0.1f) {  // Minimum drag threshold
                    dragX = (dragX / length) * MOUSE_VELOCITY_STRENGTH;
                    dragY = (dragY / length) * MOUSE_VELOCITY_STRENGTH;

                    addVelocityAtMouse(grid, currentMousePos.x, currentMousePos.y, dragX, dragY, startX, startY, cellSize);
                }

                lastMousePos = currentMousePos;
            }
        } else {
            isRightDragging = false;
        }

        // Run simulation step if running
        if (simulationRunning) {
            // Full fluid simulation (both velocity and density)
            grid.vel_step(viscosity, timeStep);
            grid.density_step(diffusionRate, timeStep);
        }

        // Render grid components
        drawCellBackgrounds(draw_list, grid, cellSize, startX, startY);
        drawGridLines(draw_list, grid, cellSize, startX, startY, gridWidth, gridHeight);

        // Draw velocity arrows if enabled
        if (showArrows) {
            drawFlowVectors(draw_list, grid, cellSize, startX, startY);
        }

        // Render UI
        renderControlPanel(grid, simulationRunning, timeStep, diffusionRate, viscosity,
                          showArrows, io, startX, startY, cellSize);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}