#ifndef GRID3D_H
#define GRID3D_H

#include <vector>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <omp.h>

class Grid3D {
public:
    int N; // Grid size (N x N x N) - number of pressure cells
    float dt;
    float diff;
    float visc; // viscosity value, high visc -> honey, low visc -> water/air
    float h = 1.0f;  // Cell size in world units
    float inv_h = 1.0f;  // 1/h

    // Pressure at cell center
    std::vector<float> p;

    // Velocities (x,y,z) (MAC/ staggered grid)
    std::vector<float> u;
    std::vector<float> v;
    std::vector<float> w;

    std::vector<float> u_old, v_old, w_old; // helper variables for vel advection

    // Density at cell center
    std::vector<float> dens;
    std::vector<float> dens_old;

    Grid3D(int size, float diffusion, float viscosity, float timestep)
        : N(size), diff(diffusion), visc(viscosity), dt(timestep)
    {
        // Pressure and density: cell centers (including ghost cells)
        int cell_count = (N + 2) * (N + 2) * (N + 2);
        p.resize(cell_count, 0.0f);
        dens.resize(cell_count, 0.0f);
        dens_old.resize(cell_count, 0.0f);

        // u velocity: at x-faces (i+1/2, j, k)
        int u_count = (N + 1) * (N + 2) * (N + 2);
        u.resize(u_count, 0.0f);
        u_old.resize(u_count, 0.0f);

        // v velocity: at y-faces (i, j+1/2, k)
        int v_count = (N + 2) * (N + 1) * (N + 2);
        v.resize(v_count, 0.0f);
        v_old.resize(v_count, 0.0f);

        // w velocity: at z-faces (i, j, k+1/2)
        int w_count = (N + 2) * (N + 2) * (N + 1);
        w.resize(w_count, 0.0f);
        w_old.resize(w_count, 0.0f);
        h = 1.0f / N;  // Cell size
        inv_h = N;      // 1/h
    }
    
    // Pressure/density at cell center (i, j, k)
    inline int P_IX(int i, int j, int k) const {
        return i + (N + 2) * j + (N + 2) * (N + 2) * k;
    }

    // u at x-face
    inline int U_IX(int i, int j, int k) const {
        return i + (N + 1) * j + (N + 1) * (N + 2) * k;
    }

    // v at y-face
    inline int V_IX(int i, int j, int k) const {
        return i + (N + 2) * j + (N + 2) * (N + 1) * k;
    }

    // w at z-face
    inline int W_IX(int i, int j, int k) const {
        return i + (N + 2) * j + (N + 2) * (N + 2) * k;
    }


    void add_density(int i, int j, int k, float amount) {
        dens[P_IX(i, j, k)] += amount;
    }

    void add_velocity(int i, int j, int k, float amountU, float amountV, float amountW) {

        // u component: average of two x-faces
        u[U_IX(i, j, k)] += amountU * 0.5f;
        u[U_IX(i+1, j, k)] += amountU * 0.5f;

        // v component: average of two y-faces
        v[V_IX(i, j, k)] += amountV * 0.5f;
        v[V_IX(i, j+1, k)] += amountV * 0.5f;

        // w component: average of two z-faces
        w[W_IX(i, j, k)] += amountW * 0.5f;
        w[W_IX(i, j, k+1)] += amountW * 0.5f;
    }

    void step() {
        vel_step();
        dens_step();
    }

    std::vector<float> getVelocityAtCellCenter(int i, int j, int k) const { // helper function that uses staggered grid for accessing each cell
        float u_avg = 0.5f * (u[U_IX(i-1, j, k)] + u[U_IX(i, j, k)]);
        float v_avg = 0.5f * (v[V_IX(i, j-1, k)] + v[V_IX(i, j, k)]);
        float w_avg = 0.5f * (w[W_IX(i, j, k-1)] + w[W_IX(i, j, k)]);

        return {u_avg, v_avg, w_avg};
    }

    void clearAllVelocities() {
        std::fill(u.begin(), u.end(), 0.0f);
        std::fill(v.begin(), v.end(), 0.0f);
        std::fill(w.begin(), w.end(), 0.0f);
        std::fill(u_old.begin(), u_old.end(), 0.0f);
        std::fill(v_old.begin(), v_old.end(), 0.0f);
        std::fill(w_old.begin(), w_old.end(), 0.0f);
    }

    void initRandomStaggeredVelocities(float magnitude = 0.5f)  {
        clearAllVelocities();

        // Initialize u at x-faces
        for (int k = 0; k <= N+1; k++) {
            for (int j = 0; j <= N+1; j++) {
                for (int i = 0; i <= N; i++) { // Note: i from 0 to N (N+1 faces)
                    u[U_IX(i, j, k)] = magnitude * ((rand() / (float)RAND_MAX) * 2.0f - 1.0f);
                }
            }
        }

        // Initialize v at y-faces
        for (int k = 0; k <= N+1; k++) {
            for (int j = 0; j <= N; j++) { // Note: j from 0 to N
                for (int i = 0; i <= N+1; i++) {
                    v[V_IX(i, j, k)] = magnitude * ((rand() / (float)RAND_MAX) * 2.0f - 1.0f);
                }
            }
        }

        // Initialize w at z-faces
        for (int k = 0; k <= N; k++) { // Note: k from 0 to N
            for (int j = 0; j <= N+1; j++) {
                for (int i = 0; i <= N+1; i++) {
                    w[W_IX(i, j, k)] = magnitude * ((rand() / (float)RAND_MAX) * 2.0f - 1.0f);
                }
            }
        }
        project(dt);
    }

    float interpolate_density(float x, float y, float z) const {
        x = std::max(0.5f, std::min((float)N + 0.5f, x));
        y = std::max(0.5f, std::min((float)N + 0.5f, y));
        z = std::max(0.5f, std::min((float)N + 0.5f, z));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);
        int k = (int)floor(z - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;
        float u = (z - 0.5f) - k;

        i = std::max(0, std::min(N+1, i));
        j = std::max(0, std::min(N+1, j));
        k = std::max(0, std::min(N+1, k));

        return (1-s)*(1-t)*(1-u)*dens_old[P_IX(i, j, k)] +
               s*(1-t)*(1-u)*dens_old[P_IX(i+1, j, k)] +
               (1-s)*t*(1-u)*dens_old[P_IX(i, j+1, k)] +
               s*t*(1-u)*dens_old[P_IX(i+1, j+1, k)] +
               (1-s)*(1-t)*u*dens_old[P_IX(i, j, k+1)] +
               s*(1-t)*u*dens_old[P_IX(i+1, j, k+1)] +
               (1-s)*t*u*dens_old[P_IX(i, j+1, k+1)] +
               s*t*u*dens_old[P_IX(i+1, j+1, k+1)];
    }

    float interpolate_u_old(float x, float y, float z) const {
    x = std::max(0.5f, std::min((float)N + 0.5f, x));
    y = std::max(0.5f, std::min((float)N + 1.5f, y));
    z = std::max(0.5f, std::min((float)N + 1.5f, z));

    int i = (int)floor(x - 0.5f);
    int j = (int)floor(y - 0.5f);
    int k = (int)floor(z - 0.5f);

    float s = (x - 0.5f) - i;
    float t = (y - 0.5f) - j;
    float uu = (z - 0.5f) - k;

    i = std::max(0, std::min(N, i));
    j = std::max(0, std::min(N+1, j));
    k = std::max(0, std::min(N+1, k));

    // Use u_old array for advection
    return (1-s)*(1-t)*(1-uu)*u_old[U_IX(i, j, k)] +
           s*(1-t)*(1-uu)*u_old[U_IX(i+1, j, k)] +
           (1-s)*t*(1-uu)*u_old[U_IX(i, j+1, k)] +
           s*t*(1-uu)*u_old[U_IX(i+1, j+1, k)] +
           (1-s)*(1-t)*uu*u_old[U_IX(i, j, k+1)] +
           s*(1-t)*uu*u_old[U_IX(i+1, j, k+1)] +
           (1-s)*t*uu*u_old[U_IX(i, j+1, k+1)] +
           s*t*uu*u_old[U_IX(i+1, j+1, k+1)];
}

    float interpolate_v_old(float x, float y, float z) const {
        x = std::max(0.5f, std::min((float)N + 1.5f, x));
        y = std::max(0.5f, std::min((float)N + 0.5f, y));
        z = std::max(0.5f, std::min((float)N + 1.5f, z));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);
        int k = (int)floor(z - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;
        float uu = (z - 0.5f) - k;

        i = std::max(0, std::min(N+1, i));
        j = std::max(0, std::min(N, j));
        k = std::max(0, std::min(N+1, k));

        // Use v_old array for advection
        return (1-s)*(1-t)*(1-uu)*v_old[V_IX(i, j, k)] +
               s*(1-t)*(1-uu)*v_old[V_IX(i+1, j, k)] +
               (1-s)*t*(1-uu)*v_old[V_IX(i, j+1, k)] +
               s*t*(1-uu)*v_old[V_IX(i+1, j+1, k)] +
               (1-s)*(1-t)*uu*v_old[V_IX(i, j, k+1)] +
               s*(1-t)*uu*v_old[V_IX(i+1, j, k+1)] +
               (1-s)*t*uu*v_old[V_IX(i, j+1, k+1)] +
               s*t*uu*v_old[V_IX(i+1, j+1, k+1)];
    }

    float interpolate_w_old(float x, float y, float z) const {
        x = std::max(0.5f, std::min((float)N + 1.5f, x));
        y = std::max(0.5f, std::min((float)N + 1.5f, y));
        z = std::max(0.5f, std::min((float)N + 0.5f, z));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);
        int k = (int)floor(z - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;
        float uu = (z - 0.5f) - k;

        i = std::max(0, std::min(N+1, i));
        j = std::max(0, std::min(N+1, j));
        k = std::max(0, std::min(N, k));

        // Use w_old array for advection
        return (1-s)*(1-t)*(1-uu)*w_old[W_IX(i, j, k)] +
               s*(1-t)*(1-uu)*w_old[W_IX(i+1, j, k)] +
               (1-s)*t*(1-uu)*w_old[W_IX(i, j+1, k)] +
               s*t*(1-uu)*w_old[W_IX(i+1, j+1, k)] +
               (1-s)*(1-t)*uu*w_old[W_IX(i, j, k+1)] +
               s*(1-t)*uu*w_old[W_IX(i+1, j, k+1)] +
               (1-s)*t*uu*w_old[W_IX(i, j+1, k+1)] +
               s*t*uu*w_old[W_IX(i+1, j+1, k+1)];
    }

private:

    void set_bnd(int b, std::vector<float>& x, int size_x, int size_y, int size_z) {
        // b indicates the type of field:
        // 0 = scalar (density, temperature)
        // 1 = u velocity
        // 2 = v velocity
        // 3 = w velocity

        // Set boundaries for each face

        // X boundaries (i = 0 and i = size_x-1)
        for (int k = 1; k < size_z-1; k++) {
            for (int j = 1; j < size_y-1; j++) {

                if (b == 1) {
                    x[0 + j*size_x + k*size_x*size_y] = 0.0f;
                } else {
                    x[0 + j*size_x + k*size_x*size_y] = x[1 + j*size_x + k*size_x*size_y];
                }
                if (b == 1) {
                    x[(size_x-1) + j*size_x + k*size_x*size_y] = 0.0f;
                } else {
                    x[(size_x-1) + j*size_x + k*size_x*size_y] = x[(size_x-2) + j*size_x + k*size_x*size_y];
                }
            }
        }

        // Y boundaries (j = 0 and j = size_y-1)
        for (int k = 1; k < size_z-1; k++) {
            for (int i = 1; i < size_x-1; i++) {
                if (b == 2) {
                    x[i + 0*size_x + k*size_x*size_y] = 0.0f;
                } else {
                    x[i + 0*size_x + k*size_x*size_y] = x[i + 1*size_x + k*size_x*size_y];
                }
                if (b == 2) {
                    x[i + (size_y-1)*size_x + k*size_x*size_y] = 0.0f;
                } else {
                    x[i + (size_y-1)*size_x + k*size_x*size_y] = x[i + (size_y-2)*size_x + k*size_x*size_y];
                }
            }
        }

        for (int j = 1; j < size_y-1; j++) {
            for (int i = 1; i < size_x-1; i++) {
                if (b == 3) {
                    x[i + j*size_x + 0*size_x*size_y] = 0.0f;
                } else {
                    x[i + j*size_x + 0*size_x*size_y] = x[i + j*size_x + 1*size_x*size_y];
                }
                if (b == 3) {
                    x[i + j*size_x + (size_z-1)*size_x*size_y] = 0.0f;
                } else {
                    x[i + j*size_x + (size_z-1)*size_x*size_y] = x[i + j*size_x + (size_z-2)*size_x*size_y];
                }
            }
        }

        // Set corners (average of adjacent faces for better stability)

        // 8 corners
        int corners[8][3] = {
            {0, 0, 0}, {size_x-1, 0, 0},
            {0, size_y-1, 0}, {size_x-1, size_y-1, 0},
            {0, 0, size_z-1}, {size_x-1, 0, size_z-1},
            {0, size_y-1, size_z-1}, {size_x-1, size_y-1, size_z-1}
        };

        for (int c = 0; c < 8; c++) {
            int i = corners[c][0];
            int j = corners[c][1];
            int k = corners[c][2];

            float sum = 0.0f;
            int count = 0;

            // Average valid neighbors
            if (i > 0) { sum += x[(i-1) + j*size_x + k*size_x*size_y]; count++; }
            if (i < size_x-1) { sum += x[(i+1) + j*size_x + k*size_x*size_y]; count++; }
            if (j > 0) { sum += x[i + (j-1)*size_x + k*size_x*size_y]; count++; }
            if (j < size_y-1) { sum += x[i + (j+1)*size_x + k*size_x*size_y]; count++; }
            if (k > 0) { sum += x[i + j*size_x + (k-1)*size_x*size_y]; count++; }
            if (k < size_z-1) { sum += x[i + j*size_x + (k+1)*size_x*size_y]; count++; }

            if (count > 0) {
                x[i + j*size_x + k*size_x*size_y] = sum / count;
            }
        }
    }

    void advect_velocity(float dt) {
        // Save current velocities
        std::copy(u.begin(), u.end(), u_old.begin());
        std::copy(v.begin(), v.end(), v_old.begin());
        std::copy(w.begin(), w.end(), w_old.begin());

        // Advect u velocities (x-component)
        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 0; i <= N; i++) {
                    // Position of this u-face
                    float x = i + 0.5f;
                    float y = j + 0.5f;
                    float z = k + 0.5f;

                    // Get velocity at this position (average of velocities)
                    float u_vel = u_old[U_IX(i, j, k)];
                    float v_vel = 0.25f * (
                        v_old[V_IX(i, j-1, k)] + v_old[V_IX(i, j, k)] +
                        v_old[V_IX(i+1, j-1, k)] + v_old[V_IX(i+1, j, k)]
                    );
                    float w_vel = 0.25f * (
                        w_old[W_IX(i, j, k-1)] + w_old[W_IX(i, j, k)] +
                        w_old[W_IX(i+1, j, k-1)] + w_old[W_IX(i+1, j, k)]
                    );

                    // Backtrace (RK2 method)
                    float half_dt = 0.5f * dt;

                    // First evaluation at current position
                    float k1_u = u_vel;
                    float k1_v = v_vel;
                    float k1_w = w_vel;

                    // Second evaluation at mid-point
                    float mid_x = x - half_dt * k1_u * inv_h;
                    float mid_y = y - half_dt * k1_v * inv_h;
                    float mid_z = z - half_dt * k1_w * inv_h;

                    // Get velocity at mid-point (interpolate from old fields)
                    float mid_u = interpolate_u_old(mid_x, mid_y, mid_z);
                    float mid_v = interpolate_v_old(mid_x, mid_y, mid_z);
                    float mid_w = interpolate_w_old(mid_x, mid_y, mid_z);

                    // Final backtraced position (RK2)
                    float src_x = x - dt * mid_u * inv_h;
                    float src_y = y - dt * mid_v * inv_h;
                    float src_z = z - dt * mid_w * inv_h;

                    // Interpolate u from old field at backtraced position
                    u[U_IX(i, j, k)] = interpolate_u_old(src_x, src_y, src_z);
                }
            }
        }

        // Advect v velocities (y-component)
        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N; k++) {
            for (int j = 0; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    float x = i + 0.5f;
                    float y = j + 0.5f;
                    float z = k + 0.5f;

                    float u_vel = 0.25f * (
                        u_old[U_IX(i-1, j, k)] + u_old[U_IX(i, j, k)] +
                        u_old[U_IX(i-1, j+1, k)] + u_old[U_IX(i, j+1, k)]
                    );
                    float v_vel = v_old[V_IX(i, j, k)];
                    float w_vel = 0.25f * (
                        w_old[W_IX(i, j, k-1)] + w_old[W_IX(i, j, k)] +
                        w_old[W_IX(i, j+1, k-1)] + w_old[W_IX(i, j+1, k)]
                    );

                    // RK2 backtrace
                    float half_dt = 0.5f * dt;
                    float k1_u = u_vel;
                    float k1_v = v_vel;
                    float k1_w = w_vel;

                    float mid_x = x - half_dt * k1_u * inv_h;
                    float mid_y = y - half_dt * k1_v * inv_h;
                    float mid_z = z - half_dt * k1_w * inv_h;

                    float mid_u = interpolate_u_old(mid_x, mid_y, mid_z);
                    float mid_v = interpolate_v_old(mid_x, mid_y, mid_z);
                    float mid_w = interpolate_w_old(mid_x, mid_y, mid_z);

                    float src_x = x - dt * mid_u * inv_h;
                    float src_y = y - dt * mid_v * inv_h;
                    float src_z = z - dt * mid_w * inv_h;

                    v[V_IX(i, j, k)] = interpolate_v_old(src_x, src_y, src_z);
                }
            }
        }

        // Advect w velocities (z-component)
        #pragma omp parallel for collapse(3)
        for (int k = 0; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    float x = i + 0.5f;
                    float y = j + 0.5f;
                    float z = k + 0.5f;

                    float u_vel = 0.25f * (
                        u_old[U_IX(i-1, j, k)] + u_old[U_IX(i, j, k)] +
                        u_old[U_IX(i-1, j, k+1)] + u_old[U_IX(i, j, k+1)]
                    );
                    float v_vel = 0.25f * (
                        v_old[V_IX(i, j-1, k)] + v_old[V_IX(i, j, k)] +
                        v_old[V_IX(i, j-1, k+1)] + v_old[V_IX(i, j, k+1)]
                    );
                    float w_vel = w_old[W_IX(i, j, k)];

                    // RK2 backtrace
                    float half_dt = 0.5f * dt;
                    float k1_u = u_vel;
                    float k1_v = v_vel;
                    float k1_w = w_vel;

                    float mid_x = x - half_dt * k1_u * inv_h;
                    float mid_y = y - half_dt * k1_v * inv_h;
                    float mid_z = z - half_dt * k1_w * inv_h;

                    float mid_u = interpolate_u_old(mid_x, mid_y, mid_z);
                    float mid_v = interpolate_v_old(mid_x, mid_y, mid_z);
                    float mid_w = interpolate_w_old(mid_x, mid_y, mid_z);

                    float src_x = x - dt * mid_u * inv_h;
                    float src_y = y - dt * mid_v * inv_h;
                    float src_z = z - dt * mid_w * inv_h;

                    w[W_IX(i, j, k)] = interpolate_w_old(src_x, src_y, src_z);
                }
            }
        }

        // Apply boundary conditions
        set_bnd(1, u, N+1, N+2, N+2);
        set_bnd(2, v, N+2, N+1, N+2);
        set_bnd(3, w, N+2, N+2, N+1);
    }

    void advect_density(float dt) {
        // Save current density to old array
        std::copy(dens.begin(), dens.end(), dens_old.begin());

        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    // Start position at cell center
                    float x = (float)i + 0.5f;
                    float y = (float)j + 0.5f;
                    float z = (float)k + 0.5f;

                    // Get velocity at cell center
                    auto vel = getVelocityAtCellCenter(i, j, k);
                    float u_vel = vel[0];
                    float v_vel = vel[1];
                    float w_vel = vel[2];

                    // Backtrack to find source position
                    float srcX = x - dt * u_vel * inv_h;
                    float srcY = y - dt * v_vel * inv_h;
                    float srcZ = z - dt * w_vel * inv_h;

                    // Clamp to grid boundaries
                    srcX = std::max(0.5f, std::min((float)N + 0.5f, srcX));
                    srcY = std::max(0.5f, std::min((float)N + 0.5f, srcY));
                    srcZ = std::max(0.5f, std::min((float)N + 0.5f, srcZ));

                    // Interpolate density from old field
                    dens[P_IX(i, j, k)] = interpolate_density(srcX, srcY, srcZ);
                }
            }
        }
        set_bnd(0, dens, N+2, N+2, N+2);
    }

    void diffuse_velocity(float dt) {
        if (visc <= 0.0f) return;

        float a = dt * visc * inv_h * inv_h;

        // Save the current velocities (which are AFTER advection)
        std::vector<float> u_rhs = u;
        std::vector<float> v_rhs = v;
        std::vector<float> w_rhs = w;

        // Work arrays for Gauss-Seidel
        std::vector<float> u_new = u;
        std::vector<float> v_new = v;
        std::vector<float> w_new = w;

        // Solve diffusion for u using Gauss-Seidel with red-black ordering
        for (int iter = 0; iter < 10; iter++) {
            // Red cells (where i+j+k is even)
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 0; i <= N; i++) {
                        if ((i + j + k) % 2 == 0) {
                            int idx = U_IX(i, j, k);
                            float sum = 0.0f;
                            int count = 0;

                            if (i > 0) { sum += u_new[U_IX(i-1, j, k)]; count++; }
                            if (i < N) { sum += u_new[U_IX(i+1, j, k)]; count++; }
                            if (j > 1) { sum += u_new[U_IX(i, j-1, k)]; count++; }
                            if (j < N) { sum += u_new[U_IX(i, j+1, k)]; count++; }
                            if (k > 1) { sum += u_new[U_IX(i, j, k-1)]; count++; }
                            if (k < N) { sum += u_new[U_IX(i, j, k+1)]; count++; }

                            if (count > 0) {
                                u_new[idx] = (u_rhs[idx] + a * sum) / (1.0f + a * count);
                            }
                        }
                    }
                }
            }

            // Black cells (where i+j+k is odd)
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 0; i <= N; i++) {
                        if ((i + j + k) % 2 == 1) {
                            int idx = U_IX(i, j, k);
                            float sum = 0.0f;
                            int count = 0;

                            if (i > 0) { sum += u_new[U_IX(i-1, j, k)]; count++; }
                            if (i < N) { sum += u_new[U_IX(i+1, j, k)]; count++; }
                            if (j > 1) { sum += u_new[U_IX(i, j-1, k)]; count++; }
                            if (j < N) { sum += u_new[U_IX(i, j+1, k)]; count++; }
                            if (k > 1) { sum += u_new[U_IX(i, j, k-1)]; count++; }
                            if (k < N) { sum += u_new[U_IX(i, j, k+1)]; count++; }

                            if (count > 0) {
                                u_new[idx] = (u_rhs[idx] + a * sum) / (1.0f + a * count);
                            }
                        }
                    }
                }
            }
            set_bnd(1, u_new, N+1, N+2, N+2);
        }
        u = std::move(u_new);

        // Solve for v with red-black ordering
        for (int iter = 0; iter < 10; iter++) {
            // Red cells
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 0; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        if ((i + j + k) % 2 == 0) {
                            int idx = V_IX(i, j, k);
                            float sum = 0.0f;
                            int count = 0;

                            if (i > 1) { sum += v_new[V_IX(i-1, j, k)]; count++; }
                            if (i < N) { sum += v_new[V_IX(i+1, j, k)]; count++; }
                            if (j > 0) { sum += v_new[V_IX(i, j-1, k)]; count++; }
                            if (j < N) { sum += v_new[V_IX(i, j+1, k)]; count++; }
                            if (k > 1) { sum += v_new[V_IX(i, j, k-1)]; count++; }
                            if (k < N) { sum += v_new[V_IX(i, j, k+1)]; count++; }

                            if (count > 0) {
                                v_new[idx] = (v_rhs[idx] + a * sum) / (1.0f + a * count);
                            }
                        }
                    }
                }
            }

            // Black cells
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 0; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        if ((i + j + k) % 2 == 1) {
                            int idx = V_IX(i, j, k);
                            float sum = 0.0f;
                            int count = 0;

                            if (i > 1) { sum += v_new[V_IX(i-1, j, k)]; count++; }
                            if (i < N) { sum += v_new[V_IX(i+1, j, k)]; count++; }
                            if (j > 0) { sum += v_new[V_IX(i, j-1, k)]; count++; }
                            if (j < N) { sum += v_new[V_IX(i, j+1, k)]; count++; }
                            if (k > 1) { sum += v_new[V_IX(i, j, k-1)]; count++; }
                            if (k < N) { sum += v_new[V_IX(i, j, k+1)]; count++; }

                            if (count > 0) {
                                v_new[idx] = (v_rhs[idx] + a * sum) / (1.0f + a * count);
                            }
                        }
                    }
                }
            }
            set_bnd(2, v_new, N+2, N+1, N+2);
        }
        v = std::move(v_new);

        // Solve for w with red-black ordering
        for (int iter = 0; iter < 10; iter++) {
            // Red cells
            #pragma omp parallel for collapse(3)
            for (int k = 0; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        if ((i + j + k) % 2 == 0) {
                            int idx = W_IX(i, j, k);
                            float sum = 0.0f;
                            int count = 0;

                            if (i > 1) { sum += w_new[W_IX(i-1, j, k)]; count++; }
                            if (i < N) { sum += w_new[W_IX(i+1, j, k)]; count++; }
                            if (j > 1) { sum += w_new[W_IX(i, j-1, k)]; count++; }
                            if (j < N) { sum += w_new[W_IX(i, j+1, k)]; count++; }
                            if (k > 0) { sum += w_new[W_IX(i, j, k-1)]; count++; }
                            if (k < N) { sum += w_new[W_IX(i, j, k+1)]; count++; }

                            if (count > 0) {
                                w_new[idx] = (w_rhs[idx] + a * sum) / (1.0f + a * count);
                            }
                        }
                    }
                }
            }

            // Black cells
            #pragma omp parallel for collapse(3)
            for (int k = 0; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        if ((i + j + k) % 2 == 1) {
                            int idx = W_IX(i, j, k);
                            float sum = 0.0f;
                            int count = 0;

                            if (i > 1) { sum += w_new[W_IX(i-1, j, k)]; count++; }
                            if (i < N) { sum += w_new[W_IX(i+1, j, k)]; count++; }
                            if (j > 1) { sum += w_new[W_IX(i, j-1, k)]; count++; }
                            if (j < N) { sum += w_new[W_IX(i, j+1, k)]; count++; }
                            if (k > 0) { sum += w_new[W_IX(i, j, k-1)]; count++; }
                            if (k < N) { sum += w_new[W_IX(i, j, k+1)]; count++; }

                            if (count > 0) {
                                w_new[idx] = (w_rhs[idx] + a * sum) / (1.0f + a * count);
                            }
                        }
                    }
                }
            }
            set_bnd(3, w_new, N+2, N+2, N+1);
        }
        w = std::move(w_new);
    }

    void diffuse_density(float dt) {
        if (diff <= 0.0f) return;

        float a = dt * diff * inv_h * inv_h;

        // Save advected density as right-hand side
        std::vector<float> dens_rhs = dens;
        std::vector<float> dens_new = dens;

        for (int iter = 0; iter < 20; iter++) {
            // Red cells
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        if ((i + j + k) % 2 == 0) {
                            int idx = P_IX(i, j, k);
                            float sum = dens_new[P_IX(i-1, j, k)] + dens_new[P_IX(i+1, j, k)] +
                                        dens_new[P_IX(i, j-1, k)] + dens_new[P_IX(i, j+1, k)] +
                                        dens_new[P_IX(i, j, k-1)] + dens_new[P_IX(i, j, k+1)];

                            dens_new[idx] = (dens_rhs[idx] + a * sum) / (1 + 6 * a);
                        }
                    }
                }
            }

            // Black cells
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        if ((i + j + k) % 2 == 1) {
                            int idx = P_IX(i, j, k);
                            float sum = dens_new[P_IX(i-1, j, k)] + dens_new[P_IX(i+1, j, k)] +
                                        dens_new[P_IX(i, j-1, k)] + dens_new[P_IX(i, j+1, k)] +
                                        dens_new[P_IX(i, j, k-1)] + dens_new[P_IX(i, j, k+1)];

                            dens_new[idx] = (dens_rhs[idx] + a * sum) / (1 + 6 * a);
                        }
                    }
                }
            }
            set_bnd(0, dens_new, N+2, N+2, N+2);
        }
        dens = std::move(dens_new);
    }

    float computeDivergence(int i, int j, int k) const {
        // equation (4.22) in paper

        float u_right = u[U_IX(i, j, k)];
        float u_left = u[U_IX(i-1, j, k)];

        float v_top = v[V_IX(i, j, k)];
        float v_bottom = v[V_IX(i, j-1, k)];

        float w_front = w[W_IX(i, j, k)];
        float w_back = w[W_IX(i, j, k-1)];

        return inv_h * ((u_right - u_left) + (v_top - v_bottom) + (w_front - w_back));
    }

    void project(float dt) {
        std::vector<float> div((N+2)*(N+2)*(N+2), 0.0f);
        std::vector<float> pressure((N+2)*(N+2)*(N+2), 0.0f);

        // 1. Compute divergence at each cell center
        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    int idx = P_IX(i, j, k);

                    float u_right = u[U_IX(i, j, k)];
                    float u_left = u[U_IX(i-1, j, k)];
                    float v_top = v[V_IX(i, j, k)];
                    float v_bottom = v[V_IX(i, j-1, k)];
                    float w_front = w[W_IX(i, j, k)];
                    float w_back = w[W_IX(i, j, k-1)];

                    div[idx] = (u_right - u_left + v_top - v_bottom + w_front - w_back) * inv_h;
                }
            }
        }

        // Set divergence boundary conditions
        set_bnd(0, div, N+2, N+2, N+2);

        float h_squared = h * h;
        int num_iterations = 10;
        float sor_factor = 1.9f;

        // Solve Poisson equation with red-black ordering
        for (int iter = 0; iter < num_iterations; iter++) {
            // Red cells
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        if ((i + j + k) % 2 == 0) {
                            int idx = P_IX(i, j, k);

                            float p_sum = pressure[P_IX(i-1, j, k)] + pressure[P_IX(i+1, j, k)] +
                                          pressure[P_IX(i, j-1, k)] + pressure[P_IX(i, j+1, k)] +
                                          pressure[P_IX(i, j, k-1)] + pressure[P_IX(i, j, k+1)];

                            float new_p = (p_sum - h_squared * div[idx]) / 6.0f;
                            pressure[idx] = pressure[idx] + sor_factor * (new_p - pressure[idx]);
                        }
                    }
                }
            }

            // Black cells
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        if ((i + j + k) % 2 == 1) {
                            int idx = P_IX(i, j, k);

                            float p_sum = pressure[P_IX(i-1, j, k)] + pressure[P_IX(i+1, j, k)] +
                                          pressure[P_IX(i, j-1, k)] + pressure[P_IX(i, j+1, k)] +
                                          pressure[P_IX(i, j, k-1)] + pressure[P_IX(i, j, k+1)];

                            float new_p = (p_sum - h_squared * div[idx]) / 6.0f;
                            pressure[idx] = pressure[idx] + sor_factor * (new_p - pressure[idx]);
                        }
                    }
                }
            }

            set_bnd(0, pressure, N+2, N+2, N+2);
        }

        // Update u velocities
        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N-1; i++) {
                    float p_right = pressure[P_IX(i+1, j, k)];
                    float p_left = pressure[P_IX(i, j, k)];
                    float pressure_grad = (p_right - p_left) * inv_h;

                    u[U_IX(i, j, k)] -= dt * pressure_grad;
                }
            }
        }

        // Update v velocities
        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N-1; j++) {
                for (int i = 1; i <= N; i++) {
                    float p_top = pressure[P_IX(i, j+1, k)];
                    float p_bottom = pressure[P_IX(i, j, k)];
                    float pressure_grad = (p_top - p_bottom) * inv_h;

                    v[V_IX(i, j, k)] -= dt * pressure_grad;
                }
            }
        }

        // Update w velocities
        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N-1; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    float p_front = pressure[P_IX(i, j, k+1)];
                    float p_back = pressure[P_IX(i, j, k)];
                    float pressure_grad = (p_front - p_back) * inv_h;

                    w[W_IX(i, j, k)] -= dt * pressure_grad;
                }
            }
        }

        // Apply boundary conditions to velocities
        set_bnd(1, u, N+1, N+2, N+2);
        set_bnd(2, v, N+2, N+1, N+2);
        set_bnd(3, w, N+2, N+2, N+1);

        // IMPORTANT: Store pressure for visualization/debugging
        p = pressure;

        // DEBUG: Print max divergence to check if projection worked
        float max_div = 0.0f;

        #pragma omp parallel for collapse(3) reduction(max:max_div)
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    max_div = std::max(max_div, std::abs(div[P_IX(i, j, k)]));
                }
            }
        }
        if (max_div > 0.001f) {
            std::cout << "Max divergence before projection: " << max_div << std::endl;
        }
    }

    void dissipate_density(float dt, float alpha = 0.1f) {
        float factor = 1.0f / (1.0f + dt * alpha);

        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    int idx = P_IX(i, j, k);
                    dens[idx] *= factor;
                }
            }
        }
        set_bnd(0, dens, N+2, N+2, N+2);
    }

    void add_forces() {
        int centerX = N / 2;
        int centerY = 2;
        int centerZ = N / 2;

        // Add circular/spiral velocity pattern
        for (int k = centerZ-2; k <= centerZ+2; k++) {
            for (int i = centerX-2; i <= centerX+2; i++) {
                float dx = (i - centerX) * 0.3f;
                float dz = (k - centerZ) * 0.3f;
                // Swirl pattern: perpendicular to radial direction
                add_velocity(i, centerY, k, -dz, 1.5f, dx);
                add_density(i, centerY, k, 150.0f);
            }
        }
    }

    void vel_step() {
        add_forces();
        std::copy(u.begin(), u.end(), u_old.begin());
        std::copy(v.begin(), v.end(), v_old.begin());
        std::copy(w.begin(), w.end(), w_old.begin());

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