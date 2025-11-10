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
const int GRID_COLS = 10;
const int GRID_ROWS = 10;

// Window settings
const int WINDOW_WIDTH = 900;
const int WINDOW_HEIGHT = 900;
const char* WINDOW_TITLE = "Fluid Grid Simulation";

// Visual settings
const float DEFAULT_CELL_SIZE = 20.0f;      // Pixels per cell
const float MIN_CELL_SIZE = 5.0f;
const float MAX_CELL_SIZE = 100.0f;

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

// Velocity color mapping
const ImU32 MAX_INFLOW_COLOR = IM_COL32(0, 0, 255, 180);      // Blue for max inflow
const ImU32 MAX_OUTFLOW_COLOR = IM_COL32(255, 0, 0, 180);     // Red for max outflow
const ImU32 NEUTRAL_COLOR = IM_COL32(128, 128, 128, 100);     // Gray for neutral
const float VELOCITY_THRESHOLD = 0.01f;                       // Minimum velocity to show color

// UI settings
const float UI_FONT_SCALE = 1.3f;
const float UI_WINDOW_ROUNDING = 8.0f;
const float UI_FRAME_ROUNDING = 4.0f;

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

void drawCellBackgrounds(ImDrawList* draw_list, const Grid& grid, float cellSize,
                        float startX, float startY) {
    for (int i = 0; i < grid.CellCountX; ++i) {
        for (int j = 0; j < grid.CellCountY; ++j) {
            float x1 = startX + i * cellSize;
            float y1 = startY + j * cellSize;
            float x2 = x1 + cellSize;
            float y2 = y1 + cellSize;

            draw_list->AddRectFilledMultiColor(
                ImVec2(x1, y1), ImVec2(x2, y2),
                BACKGROUND_TOP, BACKGROUND_TOP,
                BACKGROUND_BOTTOM, BACKGROUND_BOTTOM
            );
        }
    }
}

void drawVelocityColors(ImDrawList* draw_list, const Grid& grid, float cellSize,
                       float startX, float startY) {
    // First pass: find max divergence magnitude for normalization
    float maxDivergence = 0.0f;
    for (int i = 0; i < grid.CellCountX; ++i) {
        for (int j = 0; j < grid.CellCountY; ++j) {
            // Calculate divergence for this cell
            // Divergence = (right_flow - left_flow) + (top_flow - bottom_flow)
            float divergence = 0.0f;

            // Horizontal divergence: flow in from left, out through right
            if (i > 0) divergence += grid.cells[i-1][j].right_force;  // Flow in from left
            if (i < grid.CellCountX - 1) divergence -= grid.cells[i][j].right_force; // Flow out through right

            // Vertical divergence: flow in from bottom, out through top
            if (j > 0) divergence += grid.cells[i][j-1].top_force;    // Flow in from bottom
            if (j < grid.CellCountY - 1) divergence -= grid.cells[i][j].top_force; // Flow out through top

            maxDivergence = std::max(maxDivergence, std::abs(divergence));
        }
    }

    // Avoid division by zero
    if (maxDivergence < VELOCITY_THRESHOLD) {
        // Draw all cells as neutral to see the grid
        for (int i = 0; i < grid.CellCountX; ++i) {
            for (int j = 0; j < grid.CellCountY; ++j) {
                float x1 = startX + i * cellSize;
                float y1 = startY + j * cellSize;
                float x2 = x1 + cellSize;
                float y2 = y1 + cellSize;
                draw_list->AddRectFilled(ImVec2(x1, y1), ImVec2(x2, y2), NEUTRAL_COLOR);
            }
        }
        return;
    }

    // Second pass: draw colored rectangles for each CELL
    for (int i = 0; i < grid.CellCountX + 1; ++i) {
        for (int j = 0; j < grid.CellCountY + 1; ++j) {
            // Calculate divergence for this cell
            float divergence = 0.0f;

            // Horizontal divergence
            if (i > 0) divergence += grid.cells[i-1][j].right_force;  // Flow in from left
            if (i < grid.CellCountX - 1) divergence -= grid.cells[i][j].right_force; // Flow out through right

            // Vertical divergence
            if (j > 0) divergence += grid.cells[i][j-1].top_force;    // Flow in from bottom
            if (j < grid.CellCountY - 1) divergence -= grid.cells[i][j].top_force; // Flow out through top

            // Normalize divergence to [-1, 1] range
            float normalizedDiv = divergence / maxDivergence;
            normalizedDiv = std::clamp(normalizedDiv, -1.0f, 1.0f);

            // Calculate cell position
            float x1 = startX + i * cellSize;
            float y1 = startY + j * cellSize;
            float x2 = x1 + cellSize;
            float y2 = y1 + cellSize;

            // Choose color based on divergence
            ImU32 color;
            if (normalizedDiv > 0) {
                // Positive divergence = source (red) - fluid flowing out
                float t = normalizedDiv;
                int r = (int)(128 + t * (255 - 128));
                int g = (int)(128 * (1 - t));
                int b = (int)(128 * (1 - t));
                color = IM_COL32(r, g, b, 180);
            } else if (normalizedDiv < 0) {
                // Negative divergence = sink (blue) - fluid flowing in
                float t = -normalizedDiv;
                int r = (int)(128 * (1 - t));
                int g = (int)(128 * (1 - t));
                int b = (int)(128 + t * (255 - 128));
                color = IM_COL32(r, g, b, 180);
            } else {
                // Zero divergence (gray)
                color = NEUTRAL_COLOR;
            }

            draw_list->AddRectFilled(ImVec2(x1, y1), ImVec2(x2, y2), color);
        }
    }
}

void drawGridLines(ImDrawList* draw_list, const Grid& grid, float cellSize,
                   float startX, float startY, float gridWidth, float gridHeight) {
    // Vertical lines
    for (int i = 0; i <= grid.CellCountX + 1; ++i) {
        float x = startX + i * cellSize;
        bool isWall = (i == 0 || i == grid.CellCountX + 1);
        float thickness = isWall ? WALL_THICKNESS : INNER_LINE_THICKNESS;
        ImU32 color = isWall ? WALL_COLOR : INNER_LINE_COLOR;

        draw_list->AddLine(
            ImVec2(x, startY),
            ImVec2(x, startY + gridHeight),
            color, thickness
        );
    }

    // Horizontal lines
    for (int i = 0; i <= grid.CellCountY + 1; ++i) {
        float y = startY + i * cellSize;
        bool isWall = (i == 0 || i == grid.CellCountY + 1);
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
    for (int i = 0; i < grid.CellCountX + 1; ++i) {
        for (int j = 0; j < grid.CellCountY + 1; ++j) {
            float cx = startX + (i + 0.5f) * cellSize;
            float cy = startY + (j + 0.5f) * cellSize;

            // Combine horizontal and vertical forces into a single vector
            float fx = grid.cells[i][j].right_force;
            float fy = grid.cells[i][j].top_force;

            // Calculate magnitude
            float magnitude = std::sqrt(fx * fx + fy * fy);

            if (magnitude > FLOW_THRESHOLD) {
                // Normalize direction
                float dirX = fx / magnitude;
                float dirY = fy / magnitude;

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

                // Calculate arrow head points (perpendicular vector for wings)
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

void renderControlPanel(Grid& grid, float& cellSize, bool& showArrows, bool& showVelocityColors, const ImGuiIO& io) {
    ImGui::SetNextWindowPos(ImVec2(30, 30), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(400, 320), ImGuiCond_FirstUseEver);
    ImGui::Begin("Fluid Simulation Controls", nullptr, ImGuiWindowFlags_NoCollapse);

    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.9f, 0.9f, 0.95f, 1.0f));

    ImGui::Text("Grid Configuration: %dx%d cells", grid.CellCountX, grid.CellCountY);
    ImGui::Separator();
    ImGui::Spacing();

    ImGui::SliderFloat("Cell Display Size", &cellSize, MIN_CELL_SIZE, MAX_CELL_SIZE, "%.0f px");
    grid.cellDisplaySize = cellSize;

    ImGui::Spacing();
    ImGui::Separator();
    ImGui::Spacing();

    ImGui::Text("Visualization Modes:");
    ImGui::Checkbox("Show Flow Arrows", &showArrows);
    ImGui::Checkbox("Show Velocity Colors", &showVelocityColors);

    // Show color legend
    if (showVelocityColors) {
        ImGui::Indent();
        ImGui::TextColored(ImVec4(0.0f, 0.0f, 1.0f, 1.0f), "■ Inflow (Blue)");
        ImGui::TextColored(ImVec4(1.0f, 0.0f, 0.0f, 1.0f), "■ Outflow (Red)");
        ImGui::Unindent();
    }

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::Button("Randomize Forces", ImVec2(180, 35))) {
        grid.initializeRandomVelocities();
    }
    ImGui::SameLine();
    if (ImGui::Button("Reset Forces", ImVec2(180, 35))) {
        grid.resetVelocities();
    }

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

    // Initialize simulation
    Grid grid(GRID_COLS, GRID_ROWS, 1.0f);
    grid.cellDisplaySize = DEFAULT_CELL_SIZE;

    bool showArrows = true;
    bool showVelocityColors = false;
    float cellSize = DEFAULT_CELL_SIZE;

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
        float gridWidth = (grid.CellCountX + 1) * cellSize;
        float gridHeight = (grid.CellCountY + 1) * cellSize;
        float startX = (display_w - gridWidth) / 2.0f;
        float startY = (display_h - gridHeight) / 2.0f;

        // Render grid components
        drawCellBackgrounds(draw_list, grid, cellSize, startX, startY);

        if (showVelocityColors) {
            drawVelocityColors(draw_list, grid, cellSize, startX, startY);
        }

        drawGridLines(draw_list, grid, cellSize, startX, startY, gridWidth, gridHeight);

        if (showArrows) {
            drawFlowVectors(draw_list, grid, cellSize, startX, startY);
        }

        // Render UI
        renderControlPanel(grid, cellSize, showArrows, showVelocityColors, io);

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