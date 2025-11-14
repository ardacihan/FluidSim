#ifndef GRID_H
#define GRID_H

#include <vector>
#include "random"

class Grid {
public:
    Grid(int cellCountX, int cellCountY, float cellSize)
            : CellCountX(cellCountX), CellCountY(cellCountY), CellSize(cellSize),
              // Separate arrays for horizontal (u) and vertical (v) velocities
              u(cellCountX + 1, std::vector<float>(cellCountY + 2, 0.0f)),  // u at vertical faces
              v(cellCountX + 2, std::vector<float>(cellCountY + 1, 0.0f)),  // v at horizontal faces
              d(cellCountX + 2, std::vector<float>(cellCountY + 2, 0.0f)),  // density at cell centers
              border_color{0.0f, 0.0f, 0.0f, 1.0f},
              cellDisplaySize(20.0f)
    {
        initializeWalls();
        initializeRandomVelocities();
        initializeDensity();
    }

    int CellCountX;
    int CellCountY;
    float CellSize;

    // Staggered grid arrays:
    // u[i][j] - horizontal velocity at right face of cell (i-1, j-1)
    // v[i][j] - vertical velocity at top face of cell (i-1, j-1)
    // d[i][j] - density at center of cell (i-1, j-1)
    std::vector<std::vector<float>> u;  // Size: (CellCountX+1) x (CellCountY+2)
    std::vector<std::vector<float>> v;  // Size: (CellCountX+2) x (CellCountY+1)
    std::vector<std::vector<float>> d;  // Size: (CellCountX+2) x (CellCountY+2)

    // Visualization
    float border_color[4];
    float cellDisplaySize;

    void initializeWalls();
    void initializeRandomVelocities();
    void resetVelocities();
    void debugBorderVelocities();
    void initializeDensity();
    void clearDensity();

    void diffuse(int b, float diff, float dt);

    void set_bnd(int b, std::vector<std::vector<float>> &x);


    // Helper functions to check if a cell is a wall
    bool isWallCellU(int i, int j) const {
        return (i == 0 || i == CellCountX || j == 0 || j == CellCountY + 1);
    }
    bool isWallCellV(int i, int j) const {
        return (i == 0 || i == CellCountX + 1 || j == 0 || j == CellCountY);
    }
    bool isWallCellP(int i, int j) const {
        return (i == 0 || i == CellCountX + 1 || j == 0 || j == CellCountY + 1);
    }

    void setBorderColor(float r, float g, float b, float a = 1.0f) {
        border_color[0] = r;
        border_color[1] = g;
        border_color[2] = b;
        border_color[3] = a;
    }

    // Updated display dimensions
    float getDisplayWidth() const { return (CellCountX + 2) * cellDisplaySize; }
    float getDisplayHeight() const { return (CellCountY + 2) * cellDisplaySize; }

    // Add density to a cell
    void add_density(int i, int j, float amount) {
        // Bounds checking
        if (i >= 0 && i < static_cast<int>(d.size()) &&
            j >= 0 && j < static_cast<int>(d[0].size())) {
            d[i][j] += amount;
        }
    }

    // Get density at a cell
    float get_density(int i, int j) const {
        if (i >= 0 && i < static_cast<int>(d.size()) &&
            j >= 0 && j < static_cast<int>(d[0].size())) {
            return d[i][j];
        }
        return 0.0f;
    }

    // Set density at a cell
    void set_density(int i, int j, float value) {
        if (i >= 0 && i < static_cast<int>(d.size()) &&
            j >= 0 && j < static_cast<int>(d[0].size())) {
            d[i][j] = value;
        }
    }

    void diffuse_step(float diff, float dt) { diffuse(0, diff, dt);  }

    void advect(int b, float dt);

    void density_step(float DIFFUSION_RATE, float TIME_STEP);

    void vel_step(float visc, float dt);

    void project();

    void advect_velocity(float dt);
};

#endif // GRID_H