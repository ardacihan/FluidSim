#include "Grid.h"
#include "math.h"
#include <random>
#include <iostream>

void Grid::initializeWalls() {
    // Set wall boundary conditions for u (no-slip)
    for (int i = 0; i <= CellCountX; ++i) {
        for (int j = 0; j <= CellCountY + 1; ++j) {
            if (isWallCellU(i, j)) {
                u[i][j] = 0.0f;
            }
        }
    }

    // Set wall boundary conditions for v (no-slip)
    for (int i = 0; i <= CellCountX + 1; ++i) {
        for (int j = 0; j <= CellCountY; ++j) {
            if (isWallCellV(i, j)) {
                v[i][j] = 0.0f;
            }
        }
    }
}

void Grid::initializeRandomVelocities() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(-1.0f, 1.0f);

    // Initialize u velocities (horizontal component at vertical faces)
    for (int i = 1; i < CellCountX; ++i) {
        for (int j = 1; j <= CellCountY; ++j) {
            u[i][j] = dis(gen);
        }
    }

    // Initialize v velocities (vertical component at horizontal faces)
    for (int i = 1; i <= CellCountX; ++i) {
        for (int j = 1; j < CellCountY; ++j) {
            v[i][j] = dis(gen);
        }
    }
}

void Grid::resetVelocities() {
    // Reset u velocities
    for (int i = 0; i <= CellCountX; ++i) {
        for (int j = 0; j <= CellCountY + 1; ++j) {
            if (!isWallCellU(i, j)) {
                u[i][j] = 0.0f;
            }
        }
    }

    // Reset v velocities
    for (int i = 0; i <= CellCountX + 1; ++i) {
        for (int j = 0; j <= CellCountY; ++j) {
            if (!isWallCellV(i, j)) {
                v[i][j] = 0.0f;
            }
        }
    }
}

void Grid::debugBorderVelocities() {
    std::cout << "=== BORDER VELOCITIES ===" << std::endl;
    std::cout << "\n--- LEFT BORDER U (i=0) ---" << std::endl;
    for (int j = 0; j <= CellCountY + 1; ++j) {
        std::cout << "u[0][" << j << "] = " << u[0][j] << std::endl;
    }
}

void Grid::diffuse(int b, float diff, float dt) {
    std::vector<std::vector<float>> x0 = (b == 0) ? d : ((b == 1) ? u : v);
    std::vector<std::vector<float>>& x = (b == 0) ? d : ((b == 1) ? u : v);

    float a = dt * diff * CellCountX * CellCountY;

    int iMax = (b == 1) ? CellCountX : CellCountX + 1;
    int jMax = (b == 2) ? CellCountY : CellCountY + 1;

    // Gauss-Seidel iteration
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i < iMax; i++) {
            for (int j = 1; j < jMax; j++) {
                x[i][j] = (x0[i][j] + a * (x[i-1][j] + x[i+1][j] +
                                          x[i][j-1] + x[i][j+1])) / (1 + 4*a);
            }
        }
        set_bnd(b, x);
    }
}

void Grid::set_bnd(int b, std::vector<std::vector<float>>& x) {
    int NX = CellCountX;
    int NY = CellCountY;

    // For density field (b=0): size is (NX+2) x (NY+2)
    if (b == 0) {
        for (int i = 1; i <= NX; ++i) {
            x[i][0] = x[i][1];           // Bottom
            x[i][NY+1] = x[i][NY];       // Top
        }
        for (int j = 1; j <= NY; ++j) {
            x[0][j] = x[1][j];           // Left
            x[NX+1][j] = x[NX][j];       // Right
        }
        // Corners
        x[0][0] = 0.5f * (x[1][0] + x[0][1]);
        x[0][NY+1] = 0.5f * (x[1][NY+1] + x[0][NY]);
        x[NX+1][0] = 0.5f * (x[NX][0] + x[NX+1][1]);
        x[NX+1][NY+1] = 0.5f * (x[NX][NY+1] + x[NX+1][NY]);
    }
    // For u velocity (b=1): size is (NX+1) x (NY+2)
    else if (b == 1) {
        for (int i = 1; i < NX; ++i) {
            x[i][0] = x[i][1];           // Bottom
            x[i][NY+1] = x[i][NY];       // Top
        }
        // Left and right walls: no-slip (already zero)
        x[0][0] = 0.0f;
        x[NX][0] = 0.0f;
        x[0][NY+1] = 0.0f;
        x[NX][NY+1] = 0.0f;
    }
    // For v velocity (b=2): size is (NX+2) x (NY+1)
    else if (b == 2) {
        for (int j = 1; j < NY; ++j) {
            x[0][j] = x[1][j];           // Left
            x[NX+1][j] = x[NX][j];       // Right
        }
        // Top and bottom walls: no-slip (already zero)
        x[0][0] = 0.0f;
        x[NX+1][0] = 0.0f;
        x[0][NY] = 0.0f;
        x[NX+1][NY] = 0.0f;
    }
}

void Grid::initializeDensity() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dis(0.0f, 1.0f);

    for (int i = 1; i <= CellCountX; ++i) {
        for (int j = 1; j <= CellCountY; ++j) {
            d[i][j] = dis(gen);
        }
    }
}

void Grid::clearDensity() {
    for (int i = 0; i <= CellCountX + 1; ++i) {
        for (int j = 0; j <= CellCountY + 1; ++j) {
            d[i][j] = 0.0f;
        }
    }
}

void Grid::advect(int b, float dt) {
    std::vector<std::vector<float>>& field = (b == 0) ? d : ((b == 1) ? u : v);
    std::vector<std::vector<float>> field0 = field;

    float dt0_x = dt * CellCountX;
    float dt0_y = dt * CellCountY;

    int iMax = (b == 1) ? CellCountX : CellCountX + 1;
    int jMax = (b == 2) ? CellCountY : CellCountY + 1;

    for (int i = 1; i < iMax; i++) {
        for (int j = 1; j < jMax; j++) {
            // Interpolate velocity at current position
            float vel_u, vel_v;

            if (b == 1) {
                // For u: already at correct location
                vel_u = field0[i][j];
                vel_v = 0.25f * (v[i][j] + v[i+1][j] + v[i][j-1] + v[i+1][j-1]);
            } else if (b == 2) {
                // For v: already at correct location
                vel_u = 0.25f * (u[i][j] + u[i-1][j] + u[i][j+1] + u[i-1][j+1]);
                vel_v = field0[i][j];
            } else {
                // For density at cell center
                vel_u = 0.5f * (u[i][j] + u[i-1][j]);
                vel_v = 0.5f * (v[i][j] + v[i][j-1]);
            }

            // Trace backwards
            float x = i - dt0_x * vel_u;
            float y = j - dt0_y * vel_v;

            // Clamp to valid range
            x = std::max(0.5f, std::min(x, (float)iMax - 0.5f));
            y = std::max(0.5f, std::min(y, (float)jMax - 0.5f));

            // Bilinear interpolation
            int i0 = (int)x, i1 = i0 + 1;
            int j0 = (int)y, j1 = j0 + 1;

            float s1 = x - i0, s0 = 1 - s1;
            float t1 = y - j0, t0 = 1 - t1;

            field[i][j] = s0 * (t0 * field0[i0][j0] + t1 * field0[i0][j1]) +
                          s1 * (t0 * field0[i1][j0] + t1 * field0[i1][j1]);
        }
    }
    set_bnd(b, field);
}

void Grid::density_step(float diff, float dt) {
    diffuse(0, diff, dt);
    advect(0, dt);
}

void Grid::vel_step(float visc, float dt) {
    // Diffuse velocity
    diffuse(1, visc, dt);
    diffuse(2, visc, dt);

    // Project to make incompressible
    project();

    // Advect velocity
    advect(1, dt);
    advect(2, dt);

    // Project again
    project();
}

void Grid::project() {
    int N = CellCountX;
    float h = 1.0f / N;

    std::vector<std::vector<float>> div(N+2, std::vector<float>(N+2, 0.0f));
    std::vector<std::vector<float>> p(N+2, std::vector<float>(N+2, 0.0f));

    // FIXED: Calculate divergence using proper staggered grid indexing
    // For cell (i,j), the divergence is:
    // div = (u_right - u_left) / h + (v_top - v_bottom) / h
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            // u[i][j] is velocity at right face of cell (i,j)
            // u[i-1][j] is velocity at left face of cell (i,j)
            // v[i][j] is velocity at top face of cell (i,j)
            // v[i][j-1] is velocity at bottom face of cell (i,j)
            div[i][j] = -0.5f * h * ((u[i][j] - u[i-1][j]) + (v[i][j] - v[i][j-1]));
            p[i][j] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);

    // Solve for pressure using Gauss-Seidel
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                p[i][j] = (div[i][j] + p[i-1][j] + p[i+1][j] +
                          p[i][j-1] + p[i][j+1]) / 4;
            }
        }
        set_bnd(0, p);
    }

    // FIXED: Subtract pressure gradient at proper locations
    // For u[i][j] (velocity at face between cells i and i+1):
    // gradient = (p[i+1][j] - p[i][j]) / h
    for (int i = 1; i < N; i++) {
        for (int j = 1; j <= N; j++) {
            u[i][j] -= 0.5f * (p[i+1][j] - p[i][j]) / h;
        }
    }

    // For v[i][j] (velocity at face between cells j and j+1):
    // gradient = (p[i][j+1] - p[i][j]) / h
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j < N; j++) {
            v[i][j] -= 0.5f * (p[i][j+1] - p[i][j]) / h;
        }
    }

    set_bnd(1, u);
    set_bnd(2, v);
}