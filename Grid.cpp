#include "Grid.h"
#include "math.h"
#include <random>
#include <iostream>

float Grid::get_diverge_at(int x, int y) {
    // We represent incoming velocity as positive and outgoing velocity as negative
    return cells[y][x].bottom_force + cells[y][x].top_force + cells[y][x].right_force + cells[y][x].left_force;
}

void Grid::initializeWalls() {
    // Mark all border cells as walls
    for (int i = 0; i < CellCountX + 2; i++) {
        for (int j = 0; j < CellCountY + 2; j++) {
            // Left and right borders
            if (i == 0 || i == CellCountX + 1) {
                cells[i][j].isWall = true;
                // Set all forces to zero for wall cells
                cells[i][j].left_force = 0.0f;
                cells[i][j].right_force = 0.0f;
                cells[i][j].top_force = 0.0f;
                cells[i][j].bottom_force = 0.0f;
            }
            // Top and bottom borders
            if (j == 0 || j == CellCountY + 1) {
                cells[i][j].isWall = true;
                // Set all forces to zero for wall cells
                cells[i][j].left_force = 0.0f;
                cells[i][j].right_force = 0.0f;
                cells[i][j].top_force = 0.0f;
                cells[i][j].bottom_force = 0.0f;
            }
        }
    }
}

bool Grid::validateNeighborhoodForces() {
    const float tolerance = 1e-5f;  // Small tolerance for floating point comparison
    bool allValid = true;

    // Check horizontal neighbors (left-right relationships)
    for (int i = 0; i < CellCountX + 1; i++) {  // Updated range
        for (int j = 0; j < CellCountY + 2; j++) {
            // Skip if either cell is a wall
            if (cells[i][j].isWall || cells[i+1][j].isWall) continue;

            if (std::abs(cells[i][j].right_force + cells[i+1][j].left_force) > tolerance) {
                allValid = false;
            }
        }
    }

    // Check vertical neighbors (top-bottom relationships)
    for (int i = 0; i < CellCountX + 2; i++) {
        for (int j = 0; j < CellCountY + 1; j++) {  // Updated range
            // Skip if either cell is a wall
            if (cells[i][j].isWall || cells[i][j+1].isWall) continue;

            if (std::abs(cells[i][j].bottom_force + cells[i][j+1].top_force) > tolerance) {
                allValid = false;
            }
        }
    }

    return allValid;
}

void Grid::initializeRandomVelocities() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-1.0f, 1.0f);

    // Horizontal forces (between columns i and i+1)
    for (int i = 0; i < CellCountX + 1; i++) {  // Updated range
        for (int j = 0; j < CellCountY + 2; j++) {
            // Skip if either cell is a wall
            if (cells[i][j].isWall || cells[i+1][j].isWall) continue;

            float right_force = dis(gen);
            cells[i][j].right_force = right_force;
            cells[i+1][j].left_force = -right_force;
        }
    }

    // Vertical forces (between rows j and j+1)
    for (int i = 0; i < CellCountX + 2; i++) {
        for (int j = 0; j < CellCountY + 1; j++) {  // Updated range
            // Skip if either cell is a wall
            if (cells[i][j].isWall || cells[i][j+1].isWall) continue;

            float top_force = dis(gen);
            cells[i][j].top_force = top_force;
            cells[i][j+1].bottom_force = -top_force;
        }
    }

    debugBorderVelocities();
}

void Grid::resetVelocities() {
    for (int i = 0; i < CellCountX + 2; i++) {
        for (int j = 0; j < CellCountY + 2; j++) {
            // Only reset fluid cells, keep wall cells at zero
            if (!cells[i][j].isWall) {
                cells[i][j].right_force = 0.0f;
                cells[i][j].left_force = 0.0f;
                cells[i][j].bottom_force = 0.0f;
                cells[i][j].top_force = 0.0f;
            }
        }
    }
}

void Grid::debugBorderVelocities() {
    std::cout << "=== BORDER VELOCITIES ===" << std::endl;

    // Check left border (i=0)
    std::cout << "\n--- LEFT BORDER (i=0) ---" << std::endl;
    for (int j = 0; j < CellCountY + 2; j++) {
        std::cout << "Left border [" << 0 << "][" << j << "]: "
                  << "left=" << cells[0][j].left_force
                  << ", right=" << cells[0][j].right_force
                  << ", bottom=" << cells[0][j].bottom_force
                  << ", top=" << cells[0][j].top_force
                  << ", isWall=" << cells[0][j].isWall << std::endl;
    }

    // Check right border (i=CellCountX+1)
    std::cout << "\n--- RIGHT BORDER (i=" << CellCountX + 1 << ") ---" << std::endl;
    for (int j = 0; j < CellCountY + 2; j++) {
        std::cout << "Right border [" << CellCountX + 1 << "][" << j << "]: "
                  << "left=" << cells[CellCountX + 1][j].left_force
                  << ", right=" << cells[CellCountX + 1][j].right_force
                  << ", bottom=" << cells[CellCountX + 1][j].bottom_force
                  << ", top=" << cells[CellCountX + 1][j].top_force
                  << ", isWall=" << cells[CellCountX + 1][j].isWall << std::endl;
    }

    // Check bottom border (j=0)
    std::cout << "\n--- BOTTOM BORDER (j=0) ---" << std::endl;
    for (int i = 0; i < CellCountX + 2; i++) {
        std::cout << "Bottom border [" << i << "][0]: "
                  << "left=" << cells[i][0].left_force
                  << ", right=" << cells[i][0].right_force
                  << ", bottom=" << cells[i][0].bottom_force
                  << ", top=" << cells[i][0].top_force
                  << ", isWall=" << cells[i][0].isWall << std::endl;
    }

    // Check top border (j=CellCountY+1)
    std::cout << "\n--- TOP BORDER (j=" << CellCountY + 1 << ") ---" << std::endl;
    for (int i = 0; i < CellCountX + 2; i++) {
        std::cout << "Top border [" << i << "][" << CellCountY + 1 << "]: "
                  << "left=" << cells[i][CellCountY + 1].left_force
                  << ", right=" << cells[i][CellCountY + 1].right_force
                  << ", bottom=" << cells[i][CellCountY + 1].bottom_force
                  << ", top=" << cells[i][CellCountY + 1].top_force
                  << ", isWall=" << cells[i][CellCountY + 1].isWall << std::endl;
    }

    // Check if any wall cells have non-zero velocities
    bool hasLeakingWalls = false;
    for (int i = 0; i < CellCountX + 2; i++) {
        for (int j = 0; j < CellCountY + 2; j++) {
            if (cells[i][j].isWall &&
                (cells[i][j].left_force != 0.0f || cells[i][j].right_force != 0.0f ||
                 cells[i][j].top_force != 0.0f || cells[i][j].bottom_force != 0.0f)) {
                std::cout << "LEAK: Wall cell [" << i << "][" << j << "] has non-zero force!" << std::endl;
                hasLeakingWalls = true;
            }
        }
    }

    std::cout << "\nWall leakage: " << (hasLeakingWalls ? "YES - WALLS ARE LEAKING!" : "NO - All good!") << std::endl;
}