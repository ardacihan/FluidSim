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

int GRID_SIZE = 32;
int WINDOW_WIDTH = 1280;
int WINDOW_HEIGHT = 720;

// Fixed color scaling values
float MAX_PRESSURE = 1.0f;
float MAX_VELOCITY = 2.0f;
float DENSITY_ALPHA_SCALE = 0.8f;


// --- Global State ---
struct CameraState {
    float rotX = 30.0f;
    float rotY = -45.0f;
    float distance = 60.0f;
    ImVec2 lastMousePos;
    bool isDragging = false;
};

struct SimulationState {
    float dt = 0.01f;
    float diff = 0.001f;
    float visc = 0.001f;
    bool isRunning = false;

    // Visualization modes
    enum VisualizationMode {
        DENSITY,        // Show density (blue gradient)
        PRESSURE,        // Show pressure field (red/transparent)
        DYE
    };

    VisualizationMode visMode = DYE;

    // Layer visualization
    bool showSingleLayer = false;
    int visibleLayer = 7; // Center layer by default
    int layerAxis = 1; // 0=X, 1=Y, 2=Z

    // Optional: show velocity vectors as overlay
    bool showVelocityVectors = false;
    float vectorScale = 2.0f;
    int vectorSkip = 2;
    float minVelocityThreshold = 0.02f;

    //Add dye funtionalities
    ImVec2 lastInteractionPos;
    bool isAddingForce = false;
    bool addDye = true;
    bool addVelocity = true;
};



// --- Forward Declarations ---
namespace Graphics {
    void drawCube(float x, float y, float z, float size, float r, float g, float b, float alpha);
    void drawWireframeBox(float size);
    void drawArrow(float x1, float y1, float z1, float x2, float y2, float z2, float r, float g, float b);
    void setupPerspective(int width, int height);
    void setupCamera(const CameraState& camera);
    void renderScene(Grid3D& grid, const SimulationState& state);
    void drawCell(Grid3D& grid, const SimulationState& state, int i, int j, int k);
    void drawVelocityVectors(Grid3D& grid, const SimulationState& state);
    void drawWireframeSphere(float x, float y, float z, float radius, int segments);
    void drawSolidSphere(float x, float y, float z, float radius, int segments);
}

namespace UI {
    void setupImGui(GLFWwindow* window);
    void renderImGui(Grid3D& grid, SimulationState& state, CameraState& camera);
}

namespace Input {
    void handleCameraInput(CameraState& camera);
    void handleMouseInteraction3D(GLFWwindow* window,
                                     Grid3D& grid,
                                     SimulationState& state,
                                     const CameraState& camera);
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

void Graphics::drawWireframeSphere(float x, float y, float z, float radius, int segments = 16) {
    glPushMatrix();
    glTranslatef(x, y, z);

    glColor3f(0.8f, 0.8f, 0.2f); // Yellow sphere
    glLineWidth(2.0f);

    // Draw longitude lines
    for (int i = 0; i <= segments; i++) {
        float lat0 = M_PI * (-0.5f + (float)(i - 1) / segments);
        float z0 = radius * sin(lat0);
        float zr0 = radius * cos(lat0);

        float lat1 = M_PI * (-0.5f + (float)i / segments);
        float z1 = radius * sin(lat1);
        float zr1 = radius * cos(lat1);

        glBegin(GL_LINE_LOOP);
        for (int j = 0; j <= segments; j++) {
            float lng = 2 * M_PI * (float)(j - 1) / segments;
            float x = cos(lng);
            float y = sin(lng);

            glVertex3f(x * zr0, y * zr0, z0);
            glVertex3f(x * zr1, y * zr1, z1);
        }
        glEnd();
    }

    // Draw latitude lines
    for (int j = 0; j < segments; j++) {
        float lng0 = 2 * M_PI * (float)(j - 1) / segments;
        float lng1 = 2 * M_PI * (float)j / segments;

        glBegin(GL_LINE_LOOP);
        for (int i = 0; i <= segments; i++) {
            float lat = M_PI * (-0.5f + (float)(i - 1) / segments);
            float z = radius * sin(lat);
            float zr = radius * cos(lat);

            glVertex3f(cos(lng0) * zr, sin(lng0) * zr, z);
            glVertex3f(cos(lng1) * zr, sin(lng1) * zr, z);
        }
        glEnd();
    }

    glPopMatrix();
}

void Graphics::drawSolidSphere(float x, float y, float z, float radius, int segments = 16) {
    glPushMatrix();
    glTranslatef(x, y, z);

    // Semi-transparent yellow
    glColor4f(0.9f, 0.9f, 0.2f, 0.3f);

    // Generate sphere vertices
    std::vector<float> vertices;
    std::vector<unsigned int> indices;

    for (int i = 0; i <= segments; i++) {
        float lat = M_PI * (-0.5f + (float)i / segments);
        float z = radius * sin(lat);
        float zr = radius * cos(lat);

        for (int j = 0; j <= segments; j++) {
            float lng = 2 * M_PI * (float)j / segments;
            float x = cos(lng) * zr;
            float y = sin(lng) * zr;

            vertices.push_back(x);
            vertices.push_back(y);
            vertices.push_back(z);
        }
    }

    // Draw sphere as quads
    glBegin(GL_QUADS);
    for (int i = 0; i < segments; i++) {
        for (int j = 0; j < segments; j++) {
            int first = (i * (segments + 1)) + j;
            int second = first + segments + 1;

            glVertex3f(vertices[first*3], vertices[first*3+1], vertices[first*3+2]);
            glVertex3f(vertices[second*3], vertices[second*3+1], vertices[second*3+2]);
            glVertex3f(vertices[second*3+3], vertices[second*3+4], vertices[second*3+5]);
            glVertex3f(vertices[first*3+3], vertices[first*3+4], vertices[first*3+5]);
        }
    }
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
    glTranslatef(-(float)GRID_SIZE/2.0f, -(float)GRID_SIZE/2.0f, -(float)GRID_SIZE/2.0f);
}

void Graphics::renderScene(Grid3D& grid, const SimulationState& state) {
    // Draw Grid Boundary
    drawWireframeBox((float)GRID_SIZE);
   // drawSolidSphere(grid.sphere_x - 1.0f, grid.sphere_y - 1.0f, grid.sphere_z - 1.0f,
   //                    grid.sphere_radius, 20);
    //drawWireframeSphere(grid.sphere_x - 1.0f, grid.sphere_y - 1.0f, grid.sphere_z - 1.0f,
    //                      grid.sphere_radius, 12);

    // Draw Voxels
    if (state.showSingleLayer) {
        // Draw only one layer
        int layer = state.visibleLayer + 1; // +1 because our loops start from 1
        switch (state.layerAxis) {
            case 0: // X-axis slice
                for(int k = 1; k <= grid.N; k++) {
                    for(int j = 1; j <= grid.N; j++) {
                        int i = layer;
                        drawCell(grid, state, i, j, k);
                    }
                }
                break;
            case 1: // Y-axis slice
                for(int k = 1; k <= grid.N; k++) {
                    for(int i = 1; i <= grid.N; i++) {
                        int j = layer;
                        drawCell(grid, state, i, j, k);
                    }
                }
                break;
            case 2: // Z-axis slice
                for(int j = 1; j <= grid.N; j++) {
                    for(int i = 1; i <= grid.N; i++) {
                        int k = layer;
                        drawCell(grid, state, i, j, k);
                    }
                }
                break;
        }
    } else {
        // Draw all cells
        for(int k = 1; k <= grid.N; k++) {
            for(int j = 1; j <= grid.N; j++) {
                for(int i = 1; i <= grid.N; i++) {
                    drawCell(grid, state, i, j, k);
                }
            }
        }
    }

    // Draw velocity vectors if enabled
    if (state.showVelocityVectors) {
        drawVelocityVectors(grid, state);
    }
}

// Helper function to draw a single cell
void Graphics::drawCell(Grid3D& grid, const SimulationState& state, int i, int j, int k) {
    std::vector<float> color;
    float alpha = 0.0f;

    switch (state.visMode) {
        case SimulationState::DENSITY: {
            float d = grid.dens[grid.P_IX(i, j, k)];
            if (d > 0.05f) {
                // Blue gradient: dark blue (low) -> cyan (medium) -> white (high)
                alpha = std::min(d/100.0f, DENSITY_ALPHA_SCALE);
                float normalized = std::min(d/200.0f, 1.0f);

                if (normalized < 0.5f) {
                    // Dark blue to cyan
                    float t = normalized * 2.0f;
                    color = {0.0f, t, 1.0f, alpha};
                } else {
                    // Cyan to white
                    float t = (normalized - 0.5f) * 2.0f;
                    color = {t, 1.0f, 1.0f, alpha};
                }
            } else {
                return; // Don't draw cells with very low density
            }
            break;
        }

        case SimulationState::PRESSURE: {
            float p_val = grid.p[grid.P_IX(i, j, k)];
            float abs_p = std::abs(p_val);

            if (abs_p > 0.001f) {

                float normalized = std::min(abs_p / MAX_PRESSURE, 1.0f);


                alpha = 0.3f + 0.7f * normalized;


                float intensity = 0.5f + 0.5f * normalized;
                color = {1.0f, 0.0f, 0.0f, alpha};

            } else {
                return;
            }
            break;
        }
        case SimulationState::DYE: {
            float dye_val = grid.dye[grid.P_IX(i, j,k)];
            if (dye_val > 0.05f) {
                float normalized = std::min(dye_val / 100.0f, 1.0f);
                alpha = std::min(normalized, 0.95f);

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
            } else {
                return;
            }
            break;
        }
    }

    drawCube((float)i + 0.5f, (float)j + 0.5f, (float)k + 0.5f,
             0.85f, color[0], color[1], color[2], color[3]);
}

void Graphics::drawVelocityVectors(Grid3D& grid, const SimulationState& state) {
    if (!state.showVelocityVectors) return;

    for(int k = 1; k <= grid.N; k += state.vectorSkip) {
        for(int j = 1; j <= grid.N; j += state.vectorSkip) {
            for(int i = 1; i <= grid.N; i += state.vectorSkip) {
                if (state.showSingleLayer) {
                    int layer = state.visibleLayer + 1;
                    switch (state.layerAxis) {
                        case 0: if (i != layer) continue; break;
                        case 1: if (j != layer) continue; break;
                        case 2: if (k != layer) continue; break;
                    }
                }

                std::vector<float> vel = grid.getVelocityAtCellCenter(i, j, k);
                float vx = vel[0];
                float vy = vel[1];
                float vz = vel[2];
                float mag = std::sqrt(vx*vx + vy*vy + vz*vz);

                if (mag > state.minVelocityThreshold) {
                    float centerX = (float)i + 0.5f;
                    float centerY = (float)j + 0.5f;
                    float centerZ = (float)k + 0.5f;
                    float scale = state.vectorScale;
                    float endX = centerX + vx * scale;
                    float endY = centerY + vy * scale;
                    float endZ = centerZ + vz * scale;

                    float normalized = std::min(mag / MAX_VELOCITY, 1.0f);
                    float r = normalized;
                    float g = 0.7f;
                    float b = 1.0f - normalized;

                    drawArrow(centerX, centerY, centerZ, endX, endY, endZ, r, g, b);
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
        if (camera.distance < 5.0f) camera.distance = 5.0f;
    }
}

void Input::handleMouseInteraction3D(GLFWwindow* window,
                                     Grid3D& grid,
                                     SimulationState& state,
                                     const CameraState& camera)
{
    if (ImGui::GetIO().WantCaptureMouse) return;
    if (!state.showSingleLayer) return;

    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    bool ctrlPressed =
        glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
        glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;

    if (!ctrlPressed || !ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
        state.isAddingForce = false;
        return;
    }

    ImVec2 mousePos = ImGui::GetMousePos();

    if (!state.isAddingForce) {
        state.isAddingForce = true;
        state.lastInteractionPos = mousePos;
        return;
    }

    float dx = mousePos.x - state.lastInteractionPos.x;
    float dy = mousePos.y - state.lastInteractionPos.y;

    float velScale = 0.05f;
    float vx =  dx * velScale;
    float vy = -dy * velScale;
    float vz = 0.0f;

    float nx = mousePos.x / width;
    float ny = mousePos.y / height;

    int i = int(nx * grid.N) + 1;
    int j = int((1.0f - ny) * grid.N) + 1;
    int k = state.visibleLayer + 1;

    if (i < 1 || i > grid.N ||
        j < 1 || j > grid.N ||
        k < 1 || k > grid.N)
        return;

    switch (state.layerAxis) {
        case 0: // X slice
            std::swap(i, k);
            vx = 0.0f;
            break;
        case 1: // Y slice
            std::swap(j, k);
            vy = 0.0f;
            break;
        case 2: // Z slice
            vz = 0.0f;
            break;
    }

    if (state.addVelocity)
        grid.add_velocity(i, j, k, vx, vy, vz);

    if (state.addDye)
        grid.add_dye(i, j, k, 100.0f);

    state.lastInteractionPos = mousePos;
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

    // FPS and stats
    ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
    ImGui::Text("Grid: %d^3 cells", grid.N);

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

    // Visualization mode selection
    ImGui::Separator();
    ImGui::Text("Visualization:");

    if (ImGui::RadioButton("Density Field", state.visMode == SimulationState::DENSITY)) {
        state.visMode = SimulationState::DENSITY;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("Pressure Field", state.visMode == SimulationState::PRESSURE)) {
        state.visMode = SimulationState::PRESSURE;
    }
    ImGui::SameLine();
    if (ImGui::RadioButton("Dye", state.visMode == SimulationState::DYE)) {
        state.visMode = SimulationState::DYE;
    }

    // Show color scheme info
    ImGui::Separator();
    ImGui::Text("Color Scheme:");
    if (state.visMode == SimulationState::DENSITY) {
        ImGui::TextColored(ImVec4(0, 0, 1, 1), "Blue: Low density");
        ImGui::TextColored(ImVec4(0, 1, 1, 1), "Cyan: Medium density");
        ImGui::TextColored(ImVec4(1, 1, 1, 1), "White: High density");
    } else if(state.visMode == SimulationState::PRESSURE)  {
        ImGui::TextColored(ImVec4(1, 0, 0, 0.3f), "Transparent Red: Low pressure");
        ImGui::TextColored(ImVec4(1, 0, 0, 1.0f), "Opaque Red: High pressure");
    }else {
        ImGui::TextColored(ImVec4(0, 0, 1, 1), "Blue -> Cyan -> Green -> Yellow -> Red");
        ImGui::Text("(Based on dye concentration)");
    }

    // Vector visualization (optional overlay)
    ImGui::Separator();
    ImGui::Text("Overlay:");
    ImGui::Checkbox("Show Velocity Vectors", &state.showVelocityVectors);

    if (state.showVelocityVectors) {
        ImGui::SliderFloat("Vector Scale", &state.vectorScale, 0.1f, 20.0f);
        ImGui::SliderInt("Vector Skip", &state.vectorSkip, 1, 5);
        ImGui::SliderFloat("Min Velocity", &state.minVelocityThreshold, 0.0f, 1.0f);
    }

    // Layer visualization controls
    ImGui::Separator();
    ImGui::Text("Slice View:");
    ImGui::Checkbox("Show Single Layer", &state.showSingleLayer);

    if (state.showSingleLayer) {
        ImGui::SliderInt("Layer", &state.visibleLayer, 0, grid.N-1);

        ImGui::Text("Slice Axis:");
        ImGui::RadioButton("X", &state.layerAxis, 0);
        ImGui::SameLine();
        ImGui::RadioButton("Y", &state.layerAxis, 1);
        ImGui::SameLine();
        ImGui::RadioButton("Z", &state.layerAxis, 2);
    }

    // Velocity initialization controls
    ImGui::Separator();
    ImGui::Text("Initialization:");

    if (ImGui::Button("Clear All")) {
        grid.clearAllVelocities();
        std::fill(grid.dens.begin(), grid.dens.end(), 0.0f);
        std::fill(grid.p.begin(), grid.p.end(), 0.0f);
    }
    ImGui::SameLine();
    if (ImGui::Button("Random Velocities")) {
        grid.initRandomStaggeredVelocities(1.0f);
    }\
    ImGui::Separator();
    ImGui::Text("Slice Interaction:");
    ImGui::TextColored(ImVec4(0.5f,1.0f,0.5f,1),
        "CTRL + Left Mouse: Add to active slice");

    ImGui::Checkbox("Add Velocity", &state.addVelocity);
    ImGui::Checkbox("Add Dye", &state.addDye);


    // Camera controls
    ImGui::Separator();
    ImGui::Text("Camera:");
    ImGui::SliderFloat("Rotation X", &camera.rotX, -180.0f, 180.0f);
    ImGui::SliderFloat("Rotation Y", &camera.rotY, -180.0f, 180.0f);
    ImGui::SliderFloat("Distance", &camera.distance, 5.0f, 500.0f);

    if (ImGui::Button("Reset Camera")) {
        camera.rotX = 30.0f;
        camera.rotY = -45.0f;
        camera.distance = 15.0f;
    }

    // Statistics
    ImGui::Separator();
    ImGui::Text("Statistics:");

    // Compute some statistics
    float maxVel = 0.0f;
    float maxDens = 0.0f;
    float maxPress = 0.0f;

    for (int i = 0; i < grid.u.size(); i++) {
        if (fabs(grid.u[i]) > maxVel) maxVel = fabs(grid.u[i]);
    }
    for (int i = 0; i < grid.v.size(); i++) {
        if (fabs(grid.v[i]) > maxVel) maxVel = fabs(grid.v[i]);
    }
    for (int i = 0; i < grid.w.size(); i++) {
        if (fabs(grid.w[i]) > maxVel) maxVel = fabs(grid.w[i]);
    }

    for (int i = 0; i < grid.dens.size(); i++) {
        if (grid.dens[i] > maxDens) maxDens = grid.dens[i];
    }

    for (int i = 0; i < grid.p.size(); i++) {
        if (fabs(grid.p[i]) > maxPress) maxPress = fabs(grid.p[i]);
    }

    ImGui::Text("Max Velocity: %.4f", maxVel);
    ImGui::Text("Max Density: %.2f", maxDens);
    ImGui::Text("Max Pressure: %.4f", maxPress);

    // Reset button
    ImGui::Separator();
    if (ImGui::Button("Reset Everything")) {
        grid = Grid3D(GRID_SIZE, state.diff, state.visc, state.dt);
        state.isRunning = false;
        state.visMode = SimulationState::DENSITY;
    }

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
        WINDOW_WIDTH,
        WINDOW_HEIGHT,
        "3D Fluid Simulation - Density & Pressure",
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

    // Set up a simple light
    float lightPos[] = {10.0f, 10.0f, 10.0f, 1.0f};
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);

    float lightAmbient[] = {0.2f, 0.2f, 0.2f, 1.0f};
    float lightDiffuse[] = {0.8f, 0.8f, 0.8f, 1.0f};
    glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
}

void App::mainLoop(GLFWwindow* window, Grid3D& grid, CameraState& camera, SimulationState& state) {
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // Handle input
        Input::handleCameraInput(camera);
        Input::handleMouseInteraction3D(window, grid, state, camera);
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
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
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
int main2() {
    // Initialize window
    GLFWwindow* window = App::initializeWindow();
    if (!window) {
        return -1;
    }

    // Setup ImGui
    UI::setupImGui(window);

    // Setup graphics
    App::initializeGraphics();

    // Create simulation grid
    Grid3D grid(GRID_SIZE, 0.001f, 0.001f, 0.1f);

    // Initialize state
    CameraState camera;
    SimulationState state;

    // Run main loop
    App::mainLoop(window, grid, camera, state);

    // Cleanup
    App::shutdown(window);
    return 0;
}