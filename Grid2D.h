#ifndef GRID2D_H
#define GRID2D_H

#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <omp.h>

class Grid2D {
public:
    float sphere_radius = 6.0f;
    float sphere_x = 16.0f;
    float sphere_y = 16.0f;

    int N; // Grid size (N x N) - number of pressure cells
    float dt;
    float diff;
    float visc; // viscosity value, high visc -> honey, low visc -> water/air
    float h = 1.0f;  // Cell size in world units
    float inv_h = 1.0f;  // 1/h

    // Pressure at cell center
    std::vector<float> p;

    // Velocities (x,y) (MAC/ staggered grid)
    std::vector<float> u;
    std::vector<float> v;

    std::vector<float> u_old, v_old; // helper variables for vel advection

    // Density at cell center
    std::vector<float> dens;
    std::vector<float> dens_old;

    Grid2D(int size, float diffusion, float viscosity, float timestep)
        : N(size), diff(diffusion), visc(viscosity), dt(timestep)
    {
        // Pressure and density: cell centers (including ghost cells)
        int cell_count = (N + 2) * (N + 2);
        p.resize(cell_count, 0.0f);
        dens.resize(cell_count, 0.0f);
        dens_old.resize(cell_count, 0.0f);

        // u velocity: at x-faces (i+1/2, j)
        int u_count = (N + 1) * (N + 2);
        u.resize(u_count, 0.0f);
        u_old.resize(u_count, 0.0f);

        // v velocity: at y-faces (i, j+1/2)
        int v_count = (N + 2) * (N + 1);
        v.resize(v_count, 0.0f);
        v_old.resize(v_count, 0.0f);

        h = 1.0f / N;  // Cell size
        inv_h = N;      // 1/h
    }

    // Pressure/density at cell center (i, j)
    inline int P_IX(int i, int j) const {
        return i + (N + 2) * j;
    }

    // u at x-face
    inline int U_IX(int i, int j) const {
        return i + (N + 1) * j;
    }

    // v at y-face
    inline int V_IX(int i, int j) const {
        return i + (N + 2) * j;
    }

    void add_density(int i, int j, float amount) {
        dens[P_IX(i, j)] += amount;
    }

    void add_velocity(int i, int j, float amountU, float amountV) {
        // u component: average of two x-faces
        u[U_IX(i, j)] += amountU * 0.5f;
        u[U_IX(i+1, j)] += amountU * 0.5f;

        // v component: average of two y-faces
        v[V_IX(i, j)] += amountV * 0.5f;
        v[V_IX(i, j+1)] += amountV * 0.5f;
    }

    void step() {
        vel_step();
        dens_step();
    }

    std::vector<float> getVelocityAtCellCenter(int i, int j) const { // helper function that uses staggered grid for accessing each cell
        float u_avg = 0.5f * (u[U_IX(i-1, j)] + u[U_IX(i, j)]);
        float v_avg = 0.5f * (v[V_IX(i, j-1)] + v[V_IX(i, j)]);

        return {u_avg, v_avg};
    }

    void clearAllVelocities() {
        std::fill(u.begin(), u.end(), 0.0f);
        std::fill(v.begin(), v.end(), 0.0f);
        std::fill(u_old.begin(), u_old.end(), 0.0f);
        std::fill(v_old.begin(), v_old.end(), 0.0f);
    }

    void initRandomStaggeredVelocities(float magnitude = 0.5f)  {
        clearAllVelocities();

        // Initialize u at x-faces
        for (int j = 0; j <= N+1; j++) {
            for (int i = 0; i <= N; i++) { // Note: i from 0 to N (N+1 faces)
                u[U_IX(i, j)] = magnitude * ((rand() / (float)RAND_MAX) * 2.0f - 1.0f);
            }
        }

        // Initialize v at y-faces
        for (int j = 0; j <= N; j++) { // Note: j from 0 to N
            for (int i = 0; i <= N+1; i++) {
                v[V_IX(i, j)] = magnitude * ((rand() / (float)RAND_MAX) * 2.0f - 1.0f);
            }
        }

        project(dt);
    }

    float interpolate_density(float x, float y) const {
        x = std::max(0.5f, std::min((float)N + 0.5f, x));
        y = std::max(0.5f, std::min((float)N + 0.5f, y));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;

        i = std::max(0, std::min(N+1, i));
        j = std::max(0, std::min(N+1, j));

        return (1-s)*(1-t)*dens_old[P_IX(i, j)] +
               s*(1-t)*dens_old[P_IX(i+1, j)] +
               (1-s)*t*dens_old[P_IX(i, j+1)] +
               s*t*dens_old[P_IX(i+1, j+1)];
    }

    float interpolate_u_old(float x, float y) const {
        x = std::max(0.5f, std::min((float)N + 0.5f, x));
        y = std::max(0.5f, std::min((float)N + 1.5f, y));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;

        i = std::max(0, std::min(N, i));
        j = std::max(0, std::min(N+1, j));

        // Use u_old array for advection
        return (1-s)*(1-t)*u_old[U_IX(i, j)] +
               s*(1-t)*u_old[U_IX(i+1, j)] +
               (1-s)*t*u_old[U_IX(i, j+1)] +
               s*t*u_old[U_IX(i+1, j+1)];
    }

    float interpolate_v_old(float x, float y) const {
        x = std::max(0.5f, std::min((float)N + 1.5f, x));
        y = std::max(0.5f, std::min((float)N + 0.5f, y));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;

        i = std::max(0, std::min(N+1, i));
        j = std::max(0, std::min(N, j));

        // Use v_old array for advection
        return (1-s)*(1-t)*v_old[V_IX(i, j)] +
               s*(1-t)*v_old[V_IX(i+1, j)] +
               (1-s)*t*v_old[V_IX(i, j+1)] +
               s*t*v_old[V_IX(i+1, j+1)];
    }

private:

    void set_bnd(int b, std::vector<float>& x, int size_x, int size_y) {
        // b indicates the type of field:
        // 0 = scalar (density, temperature)
        // 1 = u velocity
        // 2 = v velocity

        // X boundaries (i = 0 and i = size_x-1)
        for (int j = 1; j < size_y-1; j++) {
            if (b == 1) {
                x[0 + j*size_x] = 0.0f;
            } else {
                x[0 + j*size_x] = x[1 + j*size_x];
            }
            if (b == 1) {
                x[(size_x-1) + j*size_x] = 0.0f;
            } else {
                x[(size_x-1) + j*size_x] = x[(size_x-2) + j*size_x];
            }
        }

        // Y boundaries (j = 0 and j = size_y-1)
        for (int i = 1; i < size_x-1; i++) {
            if (b == 2) {
                x[i + 0*size_x] = 0.0f;
            } else {
                x[i + 0*size_x] = x[i + 1*size_x];
            }
            if (b == 2) {
                x[i + (size_y-1)*size_x] = 0.0f;
            } else {
                x[i + (size_y-1)*size_x] = x[i + (size_y-2)*size_x];
            }
        }

        // Set corners (average of adjacent faces for better stability)
        int corners[4][2] = {
            {0, 0}, {size_x-1, 0},
            {0, size_y-1}, {size_x-1, size_y-1}
        };

        for (int c = 0; c < 4; c++) {
            int i = corners[c][0];
            int j = corners[c][1];

            float sum = 0.0f;
            int count = 0;

            // Average valid neighbors
            if (i > 0) { sum += x[(i-1) + j*size_x]; count++; }
            if (i < size_x-1) { sum += x[(i+1) + j*size_x]; count++; }
            if (j > 0) { sum += x[i + (j-1)*size_x]; count++; }
            if (j < size_y-1) { sum += x[i + (j+1)*size_x]; count++; }

            if (count > 0) {
                x[i + j*size_x] = sum / count;
            }
        }
    }

    void advect_velocity(float dt) {
        // Save current velocities
        std::copy(u.begin(), u.end(), u_old.begin());
        std::copy(v.begin(), v.end(), v_old.begin());

        // Advect u velocities (x-component)
        #pragma omp parallel for collapse(2)
        for (int j = 1; j <= N; j++) {
            for (int i = 0; i <= N; i++) {
                // Position of this u-face
                float x = i + 0.5f;
                float y = j + 0.5f;

                // Get velocity at this position (average of velocities)
                float u_vel = u_old[U_IX(i, j)];
                float v_vel = 0.5f * (
                    v_old[V_IX(i, j-1)] + v_old[V_IX(i, j)] +
                    v_old[V_IX(i+1, j-1)] + v_old[V_IX(i+1, j)]
                ) * 0.5f;

                // Backtrace (RK2 method)
                float half_dt = 0.5f * dt;

                // First evaluation at current position
                float k1_u = u_vel;
                float k1_v = v_vel;

                // Second evaluation at mid-point
                float mid_x = x - half_dt * k1_u * inv_h;
                float mid_y = y - half_dt * k1_v * inv_h;

                // Get velocity at mid-point (interpolate from old fields)
                float mid_u = interpolate_u_old(mid_x, mid_y);
                float mid_v = interpolate_v_old(mid_x, mid_y);

                // Final backtraced position (RK2)
                float src_x = x - dt * mid_u * inv_h;
                float src_y = y - dt * mid_v * inv_h;

                // Interpolate u from old field at backtraced position
                u[U_IX(i, j)] = interpolate_u_old(src_x, src_y);
            }
        }

        // Advect v velocities (y-component)
        #pragma omp parallel for collapse(2)
        for (int j = 0; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                float x = i + 0.5f;
                float y = j + 0.5f;

                float u_vel = 0.5f * (
                    u_old[U_IX(i-1, j)] + u_old[U_IX(i, j)] +
                    u_old[U_IX(i-1, j+1)] + u_old[U_IX(i, j+1)]
                ) * 0.5f;
                float v_vel = v_old[V_IX(i, j)];

                // RK2 backtrace
                float half_dt = 0.5f * dt;
                float k1_u = u_vel;
                float k1_v = v_vel;

                float mid_x = x - half_dt * k1_u * inv_h;
                float mid_y = y - half_dt * k1_v * inv_h;

                float mid_u = interpolate_u_old(mid_x, mid_y);
                float mid_v = interpolate_v_old(mid_x, mid_y);

                float src_x = x - dt * mid_u * inv_h;
                float src_y = y - dt * mid_v * inv_h;

                v[V_IX(i, j)] = interpolate_v_old(src_x, src_y);
            }
        }

        // Apply boundary conditions
        set_bnd(1, u, N+1, N+2);
        set_bnd(2, v, N+2, N+1);
    }

    void advect_density(float dt) {
        // Save current density to old array
        std::copy(dens.begin(), dens.end(), dens_old.begin());

        #pragma omp parallel for collapse(2)
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                // Start position at cell center
                float x = (float)i + 0.5f;
                float y = (float)j + 0.5f;

                // Get velocity at cell center
                auto vel = getVelocityAtCellCenter(i, j);
                float u_vel = vel[0];
                float v_vel = vel[1];

                // Backtrack to find source position
                float srcX = x - dt * u_vel * inv_h;
                float srcY = y - dt * v_vel * inv_h;

                // Clamp to grid boundaries
                srcX = std::max(0.5f, std::min((float)N + 0.5f, srcX));
                srcY = std::max(0.5f, std::min((float)N + 0.5f, srcY));

                // Interpolate density from old field
                dens[P_IX(i, j)] = interpolate_density(srcX, srcY);
            }
        }
        set_bnd(0, dens, N+2, N+2);
    }

    void diffuse_velocity(float dt) {
        if (visc <= 0.0f) return;

        float a = dt * visc * inv_h * inv_h;

        // Save the current velocities (which are AFTER advection)
        std::vector<float> u_rhs = u;
        std::vector<float> v_rhs = v;

        // Work arrays for Gauss-Seidel
        std::vector<float> u_new = u;
        std::vector<float> v_new = v;

        // Solve diffusion for u using Gauss-Seidel with red-black ordering
        for (int iter = 0; iter < 10; iter++) {
            // Red cells (where i+j is even)
            #pragma omp parallel for collapse(2)
            for (int j = 1; j <= N; j++) {
                for (int i = 0; i <= N; i++) {
                    if ((i + j) % 2 == 0) {
                        int idx = U_IX(i, j);
                        float sum = 0.0f;
                        int count = 0;

                        if (i > 0) { sum += u_new[U_IX(i-1, j)]; count++; }
                        if (i < N) { sum += u_new[U_IX(i+1, j)]; count++; }
                        if (j > 1) { sum += u_new[U_IX(i, j-1)]; count++; }
                        if (j < N) { sum += u_new[U_IX(i, j+1)]; count++; }

                        if (count > 0) {
                            u_new[idx] = (u_rhs[idx] + a * sum) / (1.0f + a * count);
                        }
                    }
                }
            }

            // Black cells (where i+j is odd)
            #pragma omp parallel for collapse(2)
            for (int j = 1; j <= N; j++) {
                for (int i = 0; i <= N; i++) {
                    if ((i + j) % 2 == 1) {
                        int idx = U_IX(i, j);
                        float sum = 0.0f;
                        int count = 0;

                        if (i > 0) { sum += u_new[U_IX(i-1, j)]; count++; }
                        if (i < N) { sum += u_new[U_IX(i+1, j)]; count++; }
                        if (j > 1) { sum += u_new[U_IX(i, j-1)]; count++; }
                        if (j < N) { sum += u_new[U_IX(i, j+1)]; count++; }

                        if (count > 0) {
                            u_new[idx] = (u_rhs[idx] + a * sum) / (1.0f + a * count);
                        }
                    }
                }
            }
            set_bnd(1, u_new, N+1, N+2);
        }
        u = std::move(u_new);

        // Solve for v with red-black ordering
        for (int iter = 0; iter < 10; iter++) {
            // Red cells
            #pragma omp parallel for collapse(2)
            for (int j = 0; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    if ((i + j) % 2 == 0) {
                        int idx = V_IX(i, j);
                        float sum = 0.0f;
                        int count = 0;

                        if (i > 1) { sum += v_new[V_IX(i-1, j)]; count++; }
                        if (i < N) { sum += v_new[V_IX(i+1, j)]; count++; }
                        if (j > 0) { sum += v_new[V_IX(i, j-1)]; count++; }
                        if (j < N) { sum += v_new[V_IX(i, j+1)]; count++; }

                        if (count > 0) {
                            v_new[idx] = (v_rhs[idx] + a * sum) / (1.0f + a * count);
                        }
                    }
                }
            }

            // Black cells
            #pragma omp parallel for collapse(2)
            for (int j = 0; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    if ((i + j) % 2 == 1) {
                        int idx = V_IX(i, j);
                        float sum = 0.0f;
                        int count = 0;

                        if (i > 1) { sum += v_new[V_IX(i-1, j)]; count++; }
                        if (i < N) { sum += v_new[V_IX(i+1, j)]; count++; }
                        if (j > 0) { sum += v_new[V_IX(i, j-1)]; count++; }
                        if (j < N) { sum += v_new[V_IX(i, j+1)]; count++; }

                        if (count > 0) {
                            v_new[idx] = (v_rhs[idx] + a * sum) / (1.0f + a * count);
                        }
                    }
                }
            }
            set_bnd(2, v_new, N+2, N+1);
        }
        v = std::move(v_new);
    }

    void diffuse_density(float dt) {
        if (diff <= 0.0f) return;

        float a = dt * diff * inv_h * inv_h;

        // Save advected density as right-hand side
        std::vector<float> dens_rhs = dens;
        std::vector<float> dens_new = dens;

        for (int iter = 0; iter < 20; iter++) {
            // Red cells
            #pragma omp parallel for collapse(2)
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    if ((i + j) % 2 == 0) {
                        int idx = P_IX(i, j);
                        float sum = dens_new[P_IX(i-1, j)] + dens_new[P_IX(i+1, j)] +
                                    dens_new[P_IX(i, j-1)] + dens_new[P_IX(i, j+1)];

                        dens_new[idx] = (dens_rhs[idx] + a * sum) / (1 + 4 * a);
                    }
                }
            }

            // Black cells
            #pragma omp parallel for collapse(2)
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    if ((i + j) % 2 == 1) {
                        int idx = P_IX(i, j);
                        float sum = dens_new[P_IX(i-1, j)] + dens_new[P_IX(i+1, j)] +
                                    dens_new[P_IX(i, j-1)] + dens_new[P_IX(i, j+1)];

                        dens_new[idx] = (dens_rhs[idx] + a * sum) / (1 + 4 * a);
                    }
                }
            }
            set_bnd(0, dens_new, N+2, N+2);
        }
        dens = std::move(dens_new);
    }

    float computeDivergence(int i, int j) const {
        float u_right = u[U_IX(i, j)];
        float u_left = u[U_IX(i-1, j)];

        float v_top = v[V_IX(i, j)];
        float v_bottom = v[V_IX(i, j-1)];

        return inv_h * ((u_right - u_left) + (v_top - v_bottom));
    }

    void project(float dt) {
        std::vector<float> div((N+2)*(N+2), 0.0f);
        std::vector<float> pressure((N+2)*(N+2), 0.0f);

        // 1. Compute divergence at each cell center
        #pragma omp parallel for collapse(2)
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                int idx = P_IX(i, j);

                float u_right = u[U_IX(i, j)];
                float u_left = u[U_IX(i-1, j)];
                float v_top = v[V_IX(i, j)];
                float v_bottom = v[V_IX(i, j-1)];

                div[idx] = (u_right - u_left + v_top - v_bottom) * inv_h;
            }
        }

        // Set divergence boundary conditions
        set_bnd(0, div, N+2, N+2);

        float h_squared = h * h;
        int num_iterations = 10;
        float sor_factor = 1.9f;

        // Solve Poisson equation with red-black ordering
        for (int iter = 0; iter < num_iterations; iter++) {
            // Red cells
            #pragma omp parallel for collapse(2)
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    if ((i + j) % 2 == 0) {
                        int idx = P_IX(i, j);

                        float p_sum = pressure[P_IX(i-1, j)] + pressure[P_IX(i+1, j)] +
                                      pressure[P_IX(i, j-1)] + pressure[P_IX(i, j+1)];

                        float new_p = (p_sum - h_squared * div[idx]) / 4.0f;
                        pressure[idx] = pressure[idx] + sor_factor * (new_p - pressure[idx]);
                    }
                }
            }

            // Black cells
            #pragma omp parallel for collapse(2)
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    if ((i + j) % 2 == 1) {
                        int idx = P_IX(i, j);

                        float p_sum = pressure[P_IX(i-1, j)] + pressure[P_IX(i+1, j)] +
                                      pressure[P_IX(i, j-1)] + pressure[P_IX(i, j+1)];

                        float new_p = (p_sum - h_squared * div[idx]) / 4.0f;
                        pressure[idx] = pressure[idx] + sor_factor * (new_p - pressure[idx]);
                    }
                }
            }

            set_bnd(0, pressure, N+2, N+2);
        }

        // Update u velocities
        #pragma omp parallel for collapse(2)
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N-1; i++) {
                float p_right = pressure[P_IX(i+1, j)];
                float p_left = pressure[P_IX(i, j)];
                float pressure_grad = (p_right - p_left) * inv_h;

                u[U_IX(i, j)] -= dt * pressure_grad;
            }
        }

        // Update v velocities
        #pragma omp parallel for collapse(2)
        for (int j = 1; j <= N-1; j++) {
            for (int i = 1; i <= N; i++) {
                float p_top = pressure[P_IX(i, j+1)];
                float p_bottom = pressure[P_IX(i, j)];
                float pressure_grad = (p_top - p_bottom) * inv_h;

                v[V_IX(i, j)] -= dt * pressure_grad;
            }
        }

        // Apply boundary conditions to velocities
        set_bnd(1, u, N+1, N+2);
        set_bnd(2, v, N+2, N+1);

        // IMPORTANT: Store pressure for visualization/debugging
        p = pressure;

        // DEBUG: Print max divergence to check if projection worked
        float max_div = 0.0f;

        #pragma omp parallel for collapse(2) reduction(max:max_div)
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                max_div = std::max(max_div, std::abs(div[P_IX(i, j)]));
            }
        }
        if (max_div > 0.001f) {
            std::cout << "Max divergence before projection: " << max_div << std::endl;
        }
    }

    void dissipate_density(float dt, float alpha = 0.1f) {
        float factor = 1.0f / (1.0f + dt * alpha);

        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                int idx = P_IX(i, j);
                dens[idx] *= factor;
            }
        }
        set_bnd(0, dens, N+2, N+2);
    }

    void add_forces() {
        int centerX = N / 2;
        int centerY = 2;

        // Add velocity pattern
        for (int i = centerX-2; i <= centerX+2; i++) {
            float dx = (i - centerX) * 0.3f;
            // Simple upward flow with slight horizontal variation
            add_velocity(i, centerY, -dx * 0.5f, 1.5f);
            add_density(i, centerY, 150.0f);
        }
    }

    void vel_step() {
        add_forces();
        std::copy(u.begin(), u.end(), u_old.begin());
        std::copy(v.begin(), v.end(), v_old.begin());

        diffuse_velocity(dt);
        advect_velocity(dt);
        project(dt);
    }

    void dens_step() {
        advect_density(dt);
        diffuse_density(dt);
        dissipate_density(dt, 0.32f);
    }
};

#endif