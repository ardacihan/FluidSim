#ifndef GRID_H
#define GRID_H

#include <vector>

class Cell {
public:
    float left_force;
    float right_force;
    float top_force;
    float bottom_force;
    bool isWall;  // New: indicates if this is a wall cell

    Cell(float left, float right, float top, float bottom, bool wall = false)
    : left_force(left), right_force(right), top_force(top), bottom_force(bottom), isWall(wall) {};

    // Default constructor
    Cell() : left_force(0.0f), right_force(0.0f), top_force(0.0f), bottom_force(0.0f), isWall(false) {}
};

class Grid {
public:
    Grid(int cellCountX, int cellCountY, float cellSize)
            : CellCountX(cellCountX), CellCountY(cellCountY), CellSize(cellSize),
              // Add +2 to include wall cells on both sides
              cells(cellCountX + 2, std::vector<Cell>(cellCountY + 2, Cell(0.0f, 0.0f, 0.0f, 0.0f))),
              border_color{0.0f, 0.0f, 0.0f, 1.0f},
              cellDisplaySize(20.0f)
    {
        initializeWalls();  // Initialize wall cells first
        initializeRandomVelocities();
    }

    int CellCountX;
    int CellCountY;
    float CellSize;

    // Visualization
    float border_color[4];
    float cellDisplaySize;

    float get_diverge_at(int x, int y);

    std::vector<std::vector<Cell>> cells;
    bool validateNeighborhoodForces();
    void initializeRandomVelocities();
    void initializeWalls();  // New: initialize wall cells

    void setBorderColor(float r, float g, float b, float a = 1.0f) {
        border_color[0] = r;
        border_color[1] = g;
        border_color[2] = b;
        border_color[3] = a;
    }

    void resetVelocities();

    void debugBorderVelocities();

    // Updated display dimensions to include walls
    float getDisplayWidth() const { return (CellCountX + 2) * cellDisplaySize; }
    float getDisplayHeight() const { return (CellCountY + 2) * cellDisplaySize; }

    // Helper functions to check if a cell is a wall
    bool isWallCell(int i, int j) const { return cells[i][j].isWall; }
    bool isFluidCell(int i, int j) const { return !cells[i][j].isWall; }
};

#endif // GRID_H