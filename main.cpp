#include <GL/glew.h> // Must be before GLFW
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <iostream>
#include <vector>
#include <cmath>
#include "Grid3D.h" // Your header file

// --- Configuration ---
const int GRID_SIZE = 32; // N x N x N
const int WINDOW_WIDTH = 1280;
const int WINDOW_HEIGHT = 720;

// --- Camera State ---
float camRotX = 30.0f;
float camRotY = -45.0f;
float camDist = 80.0f;
ImVec2 lastMousePos;
bool isDragging = false;

// --- Simulation State ---
float sim_dt = 0.1f;
float sim_diff = 0.0001f;
float sim_visc = 0.0001f;
bool isRunning = true;

// Helper to draw a single cube
void drawCube(float x, float y, float z, float size, float r, float g, float b, float alpha) {
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

// Draw the bounding box wireframe
void drawWireframeBox(float size) {
    glPushMatrix();
    glScalef(size, size, size);
    glColor3f(1.0f, 1.0f, 1.0f);
    glLineWidth(2.0f);

    // Simple wireframe cube
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

GLFWwindow* init_opengl() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
       return NULL;
    }

    // Set OpenGL version
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // Needed for macOS

    GLFWwindow* window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "3D Fluid Simulation", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return NULL;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable V-Sync

    // Initialize GLEW
    glewExperimental = GL_TRUE; // Needed for core profile
    if (glewInit() != GLEW_OK) {
        glfwDestroyWindow(window);
        glfwTerminate();
        return NULL;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;

    // Initialize ImGui GLFW and OpenGL3
    if (!ImGui_ImplGlfw_InitForOpenGL(window, true)) {
        std::cerr << "Failed to initialize ImGui GLFW backend" << std::endl;
        return NULL;
    }

    // Try different GLSL versions
    const char* glsl_version = "#version 130";
    if (!ImGui_ImplOpenGL3_Init(glsl_version)) {
        glsl_version = "#version 120";
        if (!ImGui_ImplOpenGL3_Init(glsl_version)) {
            std::cerr << "Failed to initialize ImGui OpenGL3 backend" << std::endl;
            return NULL;
        }
    }

    ImGui::StyleColorsDark();
    return window;
}

int main() {

    GLFWwindow* window = init_opengl();

    // Initialize 3D Grid
    std::cout << "Initializing Grid3D..." << std::endl;
    Grid3D grid(GRID_SIZE, sim_diff, sim_visc, sim_dt);

    // OpenGL Global State
    glEnable(GL_DEPTH_TEST); // Enable depth buffering
    glEnable(GL_BLEND);      // Enable transparency
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Enable lighting for better 3D effect
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);




    std::cout << "Entering main loop..." << std::endl;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // --- Input Handling (Camera) ---
        if (!ImGui::GetIO().WantCaptureMouse) {
            if (ImGui::IsMouseDown(ImGuiMouseButton_Right)) {
                ImVec2 mousePos = ImGui::GetMousePos();
                if (!isDragging) {
                    isDragging = true;
                    lastMousePos = mousePos;
                }
                float dx = mousePos.x - lastMousePos.x;
                float dy = mousePos.y - lastMousePos.y;
                camRotY += dx * 0.5f;
                camRotX += dy * 0.5f;
                lastMousePos = mousePos;
            } else {
                isDragging = false;
            }
            camDist -= ImGui::GetIO().MouseWheel * 5.0f;
            if (camDist < 10.0f) camDist = 10.0f;
        }

        // --- Simulation Step ---
        // Update grid parameters in case slider changed
        grid.dt = sim_dt;
        grid.diff = sim_diff;
        grid.visc = sim_visc;

        if (isRunning) {
            grid.step(); // Simulate steps
        }

        // --- Rendering ---
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);

        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect = (float)display_w / (float)display_h;

        // Use gluPerspective or manual perspective
        // For now, use simple manual perspective
        float nearPlane = 1.0f;
        float farPlane = 500.0f;
        float fieldOfView = 60.0f;
        float fH = tan(fieldOfView / 360.0f * 3.14159f) * nearPlane;
        float fW = fH * aspect;

        glFrustum(-fW, fW, -fH, fH, nearPlane, farPlane);

        // Setup ModelView Matrix (Camera)
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glTranslatef(0.0f, 0.0f, -camDist);
        glRotatef(camRotX, 1.0f, 0.0f, 0.0f);
        glRotatef(camRotY, 0.0f, 1.0f, 0.0f);

        // Center the grid
        glTranslatef(-(float)GRID_SIZE/2.0f, -(float)GRID_SIZE/2.0f, -(float)GRID_SIZE/2.0f);

        // Draw Grid Boundary
        drawWireframeBox((float)GRID_SIZE);


        // Draw Voxels
        for(int k=1; k<=grid.N; k++) {
            for(int j=1; j<=grid.N; j++) {
                for(int i=1; i<=grid.N; i++) {
                    float d = grid.dens[grid.IX(i, j, k)];

                    // Only draw if density is visible
                    if (d > 0.05f) {
                        // Position logic: i,j,k corresponds to cell center
                        // 0.5 offset to center coordinates
                        // Color mapping: Blue -> Cyan -> White based on density
                        float alpha = std::min(d/100.0f, 0.8f); // Normalize alpha
                        float blue = 1.0f;
                        float green = std::min(d/100.0f, 1.0f);
                        float red = std::min(d/200.0f, 1.0f);

                        // Scale slightly < 1.0 to see gaps between cells
                        drawCube((float)i + 0.5f, (float)j + 0.5f, (float)k + 0.5f,
                                0.85f, red, green, blue, alpha);
                    }
                }
            }
        }

        // --- UI Overlay ---
        ImGui::Begin("3D Fluid Controls");
        ImGui::Text("FPS: %.1f", ImGui::GetIO().Framerate);
        ImGui::Checkbox("Running", &isRunning);
        if (ImGui::Button("Reset")) {
            grid = Grid3D(GRID_SIZE, sim_diff, sim_visc, sim_dt);
        }
        ImGui::Separator();
        ImGui::SliderFloat("Viscosity", &sim_visc, 0.0f, 0.001f, "%.5f");
        ImGui::SliderFloat("Diffusion", &sim_diff, 0.0f, 0.001f, "%.5f");
        ImGui::SliderFloat("Time Step", &sim_dt, 0.0f, 0.5f);
        ImGui::Separator();
        ImGui::End();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();

    std::cout << "Clean shutdown" << std::endl;
    return 0;
}