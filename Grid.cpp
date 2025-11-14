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
    for (int i = 1; i < CellCountX; ++i) {  // Skip walls at i=0 and i=CellCountX
        for (int j = 1; j <= CellCountY; ++j) {  // Skip walls at j=0 and j=CellCountY+1
            u[i][j] = dis(gen);
        }
    }

    // Initialize v velocities (vertical component at horizontal faces)
    for (int i = 1; i <= CellCountX; ++i) {  // Skip walls at i=0 and i=CellCountX+1
        for (int j = 1; j < CellCountY; ++j) {  // Skip walls at j=0 and j=CellCountY
            v[i][j] = dis(gen);
        }
    }

    //debugBorderVelocities();
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

    // Check left border u velocities (i=0)
    std::cout << "\n--- LEFT BORDER U (i=0) ---" << std::endl;
    for (int j = 0; j <= CellCountY + 1; ++j) {
        std::cout << "u[0][" << j << "] = " << u[0][j] << std::endl;
    }

    // Check right border u velocities (i=CellCountX)
    std::cout << "\n--- RIGHT BORDER U (i=" << CellCountX << ") ---" << std::endl;
    for (int j = 0; j <= CellCountY + 1; ++j) {
        std::cout << "u[" << CellCountX << "][" << j << "] = " << u[CellCountX][j] << std::endl;
    }

    // Check bottom border v velocities (j=0)
    std::cout << "\n--- BOTTOM BORDER V (j=0) ---" << std::endl;
    for (int i = 0; i <= CellCountX + 1; ++i) {
        std::cout << "v[" << i << "][0] = " << v[i][0] << std::endl;
    }

    // Check top border v velocities (j=CellCountY)
    std::cout << "\n--- TOP BORDER V (j=" << CellCountY << ") ---" << std::endl;
    for (int i = 0; i <= CellCountX + 1; ++i) {
        std::cout << "v[" << i << "][" << CellCountY << "] = " << v[i][CellCountY] << std::endl;
    }
}

void Grid::diffuse(int b, float diff, float dt) {

    std::vector<std::vector<float>> x0 = (b == 0) ? d : ((b == 1) ? u : v);
    std::vector<std::vector<float>>& x = (b == 0) ? d : ((b == 1) ? u : v);

    float a = dt * diff * CellCountX * CellCountY;

    // Use Gauss-Seidel iteration for stable diffusion
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i <= CellCountX; i++) {
            for (int j = 1; j <= CellCountY; j++) {
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

    // Set boundaries for each side
    for (int i = 1; i <= NX; ++i) {
        // Bottom boundary (j = 0)
        if (b == 2) { // For vertical velocity component
            x[i][0] = -x[i][1]; // No-slip condition
        } else {
            x[i][0] = x[i][1]; // Continuity for other fields
        }

        // Top boundary (j = NY + 1)
        if (b == 2) { // For vertical velocity component
            x[i][NY+1] = -x[i][NY]; // No-slip condition
        } else {
            x[i][NY+1] = x[i][NY]; // Continuity for other fields
        }
    }

    for (int j = 1; j <= NY; ++j) {
        // Left boundary (i = 0)
        if (b == 1) { // For horizontal velocity component
            x[0][j] = -x[1][j]; // No-slip condition
        } else {
            x[0][j] = x[1][j]; // Continuity for other fields
        }

        // Right boundary (i = NX + 1)
        if (b == 1) { // For horizontal velocity component
            x[NX+1][j] = -x[NX][j]; // No-slip condition
        } else {
            x[NX+1][j] = x[NX][j]; // Continuity for other fields
        }
    }

    // Handle corners by averaging
    x[0][0] = 0.5f * (x[1][0] + x[0][1]);
    x[0][NY+1] = 0.5f * (x[1][NY+1] + x[0][NY]);
    x[NX+1][0] = 0.5f * (x[NX][0] + x[NX+1][1]);
    x[NX+1][NY+1] = 0.5f * (x[NX][NY+1] + x[NX+1][NY]);
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
    std::vector<std::vector<float>>& d = (b == 0) ? this->d : ((b == 1) ? u : v);
    std::vector<std::vector<float>> d0 = d;
    std::vector<std::vector<float>>& u = this->u;
    std::vector<std::vector<float>>& v = this->v;

    float dt0_x = dt * CellCountX;
    float dt0_y = dt * CellCountY;

    for (int i = 1; i <= CellCountX; i++) {
        for (int j = 1; j <= CellCountY; j++) {
            // Trace backwards in velocity field
            float x = i - dt0_x * u[i][j];
            float y = j - dt0_y * v[i][j];

            // Clamp to grid boundaries
            x = std::max(0.5f, std::min(x, CellCountX + 0.5f));
            y = std::max(0.5f, std::min(y, CellCountY + 0.5f));

            // Integer coordinates
            int i0 = (int)x, i1 = i0 + 1;
            int j0 = (int)y, j1 = j0 + 1;

            // Interpolation weights
            float s1 = x - i0, s0 = 1 - s1;
            float t1 = y - j0, t0 = 1 - t1;

            // Bilinear interpolation
            d[i][j] = s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1]) +
                      s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);
        }
    }
    set_bnd(b, d);
}

void Grid::density_step(float diff, float dt) {
    diffuse(0, diff, dt);
    advect(0, dt);
}

void Grid::vel_step(float visc, float dt) {
    // Add source terms first (handled by mouse interaction)

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

    // Calculate divergence
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            div[i][j] = -0.5f * h * (u[i][j] - u[i-1][j] +
                                    v[i][j] - v[i][j-1]);
            p[i][j] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);

    // Solve pressure using Gauss-Seidel
    for (int k = 0; k < 20; k++) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                p[i][j] = (div[i][j] + p[i-1][j] + p[i+1][j] +
                          p[i][j-1] + p[i][j+1]) / 4;
            }
        }
        set_bnd(0, p);
    }

    // Subtract pressure gradient
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            u[i][j] -= 0.5f * (p[i+1][j] - p[i-1][j]) / h;
            v[i][j] -= 0.5f * (p[i][j+1] - p[i][j-1]) / h;
        }
    }
    set_bnd(1, u);
    set_bnd(2, v);
}

void Grid::advect_velocity(float dt) {
    int NX = CellCountX;
    int NY = CellCountY;
    float dt0_x = dt * NX;
    float dt0_y = dt * NY;

    // Advect u component (needs interpolation since u is at different locations)
    std::vector<std::vector<float>> u_old = u;
    for (int i = 1; i <= NX; ++i) {
        for (int j = 1; j <= NY; ++j) {
            // u is stored at (i, j+0.5) - need to interpolate velocities to this location
            float vel_u = u_old[i][j];
            float vel_v = 0.5f * (v[i][j] + v[i-1][j]); // Interpolate v to u location

            float x = i - dt0_x * vel_u;
            float y = (j + 0.5f) - dt0_y * vel_v; // Account for u's y-position

            if (x < 0.5f) x = 0.5f; if (x > NX + 0.5f) x = NX + 0.5f;
            if (y < 0.5f) y = 0.5f; if (y > NY + 0.5f) y = NY + 0.5f;

            int i0 = (int)x, i1 = i0 + 1;
            int j0 = (int)y, j1 = j0 + 1;

            float s1 = x - i0, s0 = 1.0f - s1;
            float t1 = y - j0, t0 = 1.0f - t1;

            u[i][j] = s0 * (t0 * u_old[i0][j0] + t1 * u_old[i0][j1]) +
                      s1 * (t0 * u_old[i1][j0] + t1 * u_old[i1][j1]);
        }
    }
    set_bnd(1, u);

    // Advect v component (needs interpolation since v is at different locations)
    std::vector<std::vector<float>> v_old = v;
    for (int i = 1; i <= NX; ++i) {
        for (int j = 1; j <= NY; ++j) {
            // v is stored at (i+0.5, j) - need to interpolate velocities to this location
            float vel_u = 0.5f * (u[i][j] + u[i][j+1]); // Interpolate u to v location
            float vel_v = v_old[i][j];

            float x = (i + 0.5f) - dt0_x * vel_u; // Account for v's x-position
            float y = j - dt0_y * vel_v;

            if (x < 0.5f) x = 0.5f; if (x > NX + 0.5f) x = NX + 0.5f;
            if (y < 0.5f) y = 0.5f; if (y > NY + 0.5f) y = NY + 0.5f;

            int i0 = (int)x, i1 = i0 + 1;
            int j0 = (int)y, j1 = j0 + 1;

            float s1 = x - i0, s0 = 1.0f - s1;
            float t1 = y - j0, t0 = 1.0f - t1;

            v[i][j] = s0 * (t0 * v_old[i0][j0] + t1 * v_old[i0][j1]) +
                      s1 * (t0 * v_old[i1][j0] + t1 * v_old[i1][j1]);
        }
    }
    set_bnd(2, v);
}