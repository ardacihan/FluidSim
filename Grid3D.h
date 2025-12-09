#ifndef GRID3D_H
#define GRID3D_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

class Grid3D {
public:
    int N; // Grid size (N x N x N) - number of pressure cells
    float dt;
    float diff;
    float visc; // viscosity value, high visc -> honey, low visc -> water/air

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

    float divergenceAtCell(int i, int j, int k) const {
        float du_dx = u[U_IX(i, j, k)] - u[U_IX(i-1, j, k)];  // Forward difference
        float dv_dy = v[V_IX(i, j, k)] - v[V_IX(i, j-1, k)];  // Forward difference
        float dw_dz = w[W_IX(i, j, k)] - w[W_IX(i, j, k-1)];  // Forward difference

        return du_dx + dv_dy + dw_dz; // since we use staggered grid we look at 2 faces per dimension
    }

    bool checkDivergence(float tolerance = 1e-3f) const {
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    float div = divergenceAtCell(i, j, k);
                    if (std::abs(div) > tolerance) {
                        return false;
                    }
                }
            }
        }
        return true;
    }

    std::vector<float> getVelocityAtCellCenter(int i, int j, int k) const { // helper function that uses staggered grid for accessing each cell
        float u_avg = 0.5f * (u[U_IX(i-1, j, k)] + u[U_IX(i, j, k)]);
        float v_avg = 0.5f * (v[V_IX(i, j-1, k)] + v[V_IX(i, j, k)]);
        float w_avg = 0.5f * (w[W_IX(i, j, k-1)] + w[W_IX(i, j, k)]);

        return {u_avg, v_avg, w_avg};
    }

    float getVelocityMagnitude(int i, int j, int k) const {
        auto vel = getVelocityAtCellCenter(i, j, k);
        return std::sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
    }

    std::vector<float> getVelocityAtUFace(int i, int j, int k) const {
        float u_val = u[U_IX(i, j, k)];

        float v_avg = 0.25f * (
            v[V_IX(i, j-1, k)] + v[V_IX(i, j, k)] +
            v[V_IX(i+1, j-1, k)] + v[V_IX(i+1, j, k)]
        );

        float w_avg = 0.25f * (
            w[W_IX(i, j, k-1)] + w[W_IX(i, j, k)] +
            w[W_IX(i+1, j, k-1)] + w[W_IX(i+1, j, k)]
        );

        return {u_val, v_avg, w_avg};
    }

    std::vector<float> getVelocityAtVFace(int i, int j, int k) const {
        float u_avg = 0.25f * (
            u[U_IX(i-1, j, k)] + u[U_IX(i, j, k)] +
            u[U_IX(i-1, j+1, k)] + u[U_IX(i, j+1, k)]
        );

        float v_val = v[V_IX(i, j, k)];

        float w_avg = 0.25f * (
            w[W_IX(i, j, k-1)] + w[W_IX(i, j, k)] +
            w[W_IX(i, j+1, k-1)] + w[W_IX(i, j+1, k)]
        );

        return {u_avg, v_val, w_avg};
    }

    std::vector<float> getVelocityAtWFace(int i, int j, int k) const {
        float u_avg = 0.25f * (
            u[U_IX(i-1, j, k)] + u[U_IX(i, j, k)] +
            u[U_IX(i-1, j, k+1)] + u[U_IX(i, j, k+1)]
        );

        float v_avg = 0.25f * (
            v[V_IX(i, j-1, k)] + v[V_IX(i, j, k)] +
            v[V_IX(i, j-1, k+1)] + v[V_IX(i, j, k+1)]
        );

        float w_val = w[W_IX(i, j, k)];

        return {u_avg, v_avg, w_val};
    }

    std::vector<float> getDivergenceColor(int i, int j, int k, float threshold = 0.01f) const {
        float div = divergenceAtCell(i, j, k);
        float absDiv = std::abs(div);

        std::vector<float> color = {0.0f, 0.0f, 0.0f, 0.0f};

        if (absDiv > threshold) {
            float intensity = std::min(absDiv * 10.0f, 1.0f);
            if (div > 0) {
                color = {0.0f, 0.0f, 1.0f, intensity}; // Blue for positive
            } else {
                color = {1.0f, 0.0f, 0.0f, intensity}; // Red for negative
            }
        }

        return color;
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
        project();
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

    float interpolate_pressure(float x, float y, float z) const {
        x = std::max(0.5f, std::min((float)N + 0.5f, x));
        y = std::max(0.5f, std::min((float)N + 0.5f, y));
        z = std::max(0.5f, std::min((float)N + 0.5f, z));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);
        int k = (int)floor(z - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;
        float uu = (z - 0.5f) - k;

        i = std::max(0, std::min(N+1, i));
        j = std::max(0, std::min(N+1, j));
        k = std::max(0, std::min(N+1, k));

        return (1-s)*(1-t)*(1-uu)*p[P_IX(i, j, k)] +
               s*(1-t)*(1-uu)*p[P_IX(i+1, j, k)] +
               (1-s)*t*(1-uu)*p[P_IX(i, j+1, k)] +
               s*t*(1-uu)*p[P_IX(i+1, j+1, k)] +
               (1-s)*(1-t)*uu*p[P_IX(i, j, k+1)] +
               s*(1-t)*uu*p[P_IX(i+1, j, k+1)] +
               (1-s)*t*uu*p[P_IX(i, j+1, k+1)] +
               s*t*uu*p[P_IX(i+1, j+1, k+1)];
    }

    float interpolate_u(float x, float y, float z) const {
        x = std::max(0.5f, std::min((float)N + 0.5f, x));
        y = std::max(0.5f, std::min((float)N + 1.5f, y));
        z = std::max(0.5f, std::min((float)N + 1.5f, z));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);
        int k = (int)floor(z - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;
        float uu = (z - 0.5f) - k;  // Renamed from u to avoid name conflict

        i = std::max(0, std::min(N, i));
        j = std::max(0, std::min(N+1, j));
        k = std::max(0, std::min(N+1, k));

        // Use current u array instead of u_old
        return (1-s)*(1-t)*(1-uu)*u[U_IX(i, j, k)] +
               s*(1-t)*(1-uu)*u[U_IX(i+1, j, k)] +
               (1-s)*t*(1-uu)*u[U_IX(i, j+1, k)] +
               s*t*(1-uu)*u[U_IX(i+1, j+1, k)] +
               (1-s)*(1-t)*uu*u[U_IX(i, j, k+1)] +
               s*(1-t)*uu*u[U_IX(i+1, j, k+1)] +
               (1-s)*t*uu*u[U_IX(i, j+1, k+1)] +
               s*t*uu*u[U_IX(i+1, j+1, k+1)];
    }

    float interpolate_v(float x, float y, float z) const {
        x = std::max(0.5f, std::min((float)N + 1.5f, x));
        y = std::max(0.5f, std::min((float)N + 0.5f, y));
        z = std::max(0.5f, std::min((float)N + 1.5f, z));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);
        int k = (int)floor(z - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;
        float uu = (z - 0.5f) - k;  // Renamed from u to avoid name conflict

        i = std::max(0, std::min(N+1, i));
        j = std::max(0, std::min(N, j));
        k = std::max(0, std::min(N+1, k));

        // Use current v array instead of v_old
        return (1-s)*(1-t)*(1-uu)*v[V_IX(i, j, k)] +
               s*(1-t)*(1-uu)*v[V_IX(i+1, j, k)] +
               (1-s)*t*(1-uu)*v[V_IX(i, j+1, k)] +
               s*t*(1-uu)*v[V_IX(i+1, j+1, k)] +
               (1-s)*(1-t)*uu*v[V_IX(i, j, k+1)] +
               s*(1-t)*uu*v[V_IX(i+1, j, k+1)] +
               (1-s)*t*uu*v[V_IX(i, j+1, k+1)] +
               s*t*uu*v[V_IX(i+1, j+1, k+1)];
    }

    float interpolate_w(float x, float y, float z) const {
        x = std::max(0.5f, std::min((float)N + 1.5f, x));
        y = std::max(0.5f, std::min((float)N + 1.5f, y));
        z = std::max(0.5f, std::min((float)N + 0.5f, z));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);
        int k = (int)floor(z - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;
        float uu = (z - 0.5f) - k;  // Renamed from u to avoid name conflict

        i = std::max(0, std::min(N+1, i));
        j = std::max(0, std::min(N+1, j));
        k = std::max(0, std::min(N, k));

        // Use current w array instead of w_old
        return (1-s)*(1-t)*(1-uu)*w[W_IX(i, j, k)] +
               s*(1-t)*(1-uu)*w[W_IX(i+1, j, k)] +
               (1-s)*t*(1-uu)*w[W_IX(i, j+1, k)] +
               s*t*(1-uu)*w[W_IX(i+1, j+1, k)] +
               (1-s)*(1-t)*uu*w[W_IX(i, j, k+1)] +
               s*(1-t)*uu*w[W_IX(i+1, j, k+1)] +
               (1-s)*t*uu*w[W_IX(i, j+1, k+1)] +
               s*t*uu*w[W_IX(i+1, j+1, k+1)];
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

    void debugVelocityDamping() {
        static int frame = 0;
        frame++;

        // Track center cell velocities
        int i = N/2, j = N/2, k = N/2;
        float u_val = 0.5f * (u[U_IX(i-1, j, k)] + u[U_IX(i, j, k)]);
        float v_val = 0.5f * (v[V_IX(i, j-1, k)] + v[V_IX(i, j, k)]);
        float w_val = 0.5f * (w[W_IX(i, j, k-1)] + w[W_IX(i, j, k)]);
        float mag = std::sqrt(u_val*u_val + v_val*v_val + w_val*w_val);

        std::cout << "Frame " << frame << ": Center velocity magnitude = " << mag << std::endl;

        // Also check total energy
        float total_u = 0, total_v = 0, total_w = 0;
        for (int idx = 0; idx < u.size(); idx++) total_u += fabs(u[idx]);
        for (int idx = 0; idx < v.size(); idx++) total_v += fabs(v[idx]);
        for (int idx = 0; idx < w.size(); idx++) total_w += fabs(w[idx]);

        std::cout << "Total |u|=" << total_u << ", |v|=" << total_v << ", |w|=" << total_w << std::endl;
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
    // Save current velocities to old arrays
    std::copy(u.begin(), u.end(), u_old.begin());
    std::copy(v.begin(), v.end(), v_old.begin());
    std::copy(w.begin(), w.end(), w_old.begin());

    // Advect u (x-faces)
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 0; i <= N; i++) {
                float x = (float)i + 0.5f;
                float y = (float)j + 0.5f;
                float z = (float)k + 0.5f;

                // Get velocity at this u-face from OLD arrays
                float u_vel = u_old[U_IX(i, j, k)];
                float v_vel = 0.25f * (
                    v_old[V_IX(i, j-1, k)] + v_old[V_IX(i, j, k)] +
                    v_old[V_IX(i+1, j-1, k)] + v_old[V_IX(i+1, j, k)]
                );
                float w_vel = 0.25f * (
                    w_old[W_IX(i, j, k-1)] + w_old[W_IX(i, j, k)] +
                    w_old[W_IX(i+1, j, k-1)] + w_old[W_IX(i+1, j, k)]
                );

                // Backtrace
                float srcX = x - dt  * u_vel;
                float srcY = y - dt  * v_vel;
                float srcZ = z - dt  * w_vel;

                // Interpolate from OLD u field
                u[U_IX(i, j, k)] = interpolate_u_old(srcX, srcY, srcZ);
            }
        }
    }
    set_bnd(1, u, N+1, N+2, N+2);

    // Advect v (y-faces)
    for (int k = 1; k <= N; k++) {
        for (int j = 0; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                float x = (float)i + 0.5f;
                float y = (float)j + 0.5f;
                float z = (float)k + 0.5f;

                float u_vel = 0.25f * (
                    u_old[U_IX(i-1, j, k)] + u_old[U_IX(i, j, k)] +
                    u_old[U_IX(i-1, j+1, k)] + u_old[U_IX(i, j+1, k)]
                );
                float v_vel = v_old[V_IX(i, j, k)];
                float w_vel = 0.25f * (
                    w_old[W_IX(i, j, k-1)] + w_old[W_IX(i, j, k)] +
                    w_old[W_IX(i, j+1, k-1)] + w_old[W_IX(i, j+1, k)]
                );

                float srcX = x - dt * N * u_vel;
                float srcY = y - dt * N * v_vel;
                float srcZ = z - dt * N * w_vel;

                v[V_IX(i, j, k)] = interpolate_v_old(srcX, srcY, srcZ);
            }
        }
    }
    set_bnd(2, v, N+2, N+1, N+2);

    // Advect w (z-faces)
    for (int k = 0; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                float x = (float)i + 0.5f;
                float y = (float)j + 0.5f;
                float z = (float)k + 0.5f;

                float u_vel = 0.25f * (
                    u_old[U_IX(i-1, j, k)] + u_old[U_IX(i, j, k)] +
                    u_old[U_IX(i-1, j, k+1)] + u_old[U_IX(i, j, k+1)]
                );
                float v_vel = 0.25f * (
                    v_old[V_IX(i, j-1, k)] + v_old[V_IX(i, j, k)] +
                    v_old[V_IX(i, j-1, k+1)] + v_old[V_IX(i, j, k+1)]
                );
                float w_vel = w_old[W_IX(i, j, k)];

                float srcX = x - dt * N * u_vel;
                float srcY = y - dt * N * v_vel;
                float srcZ = z - dt * N * w_vel;

                w[W_IX(i, j, k)] = interpolate_w_old(srcX, srcY, srcZ);
            }
        }
    }
    set_bnd(3, w, N+2, N+2, N+1);
}

    void advect_density(float dt) {
        // Save current density to old array
        std::copy(dens.begin(), dens.end(), dens_old.begin());

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
                    float srcX = x - dt * N * u_vel;  // Multiply by N to convert to grid units
                    float srcY = y - dt * N * v_vel;
                    float srcZ = z - dt * N * w_vel;

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

    float a = dt * visc * N * N;

    // FIRST: Save the current velocities (which are AFTER advection)
    // These become the right-hand side of the linear system
    std::vector<float> u_rhs = u;  // w2 in paper
    std::vector<float> v_rhs = v;
    std::vector<float> w_rhs = w;

    // Work arrays for Gauss-Seidel
    std::vector<float> u_new = u;
    std::vector<float> v_new = v;
    std::vector<float> w_new = w;

    // Solve diffusion for u using Gauss-Seidel
    for (int iter = 0; iter < 20; iter++) {
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 0; i <= N; i++) {
                    int idx = U_IX(i, j, k);
                    float sum = 0.0f;
                    int count = 0;

                    // Use u_new (current iteration's values) for neighbors
                    if (i > 0) { sum += u_new[U_IX(i-1, j, k)]; count++; }
                    if (i < N) { sum += u_new[U_IX(i+1, j, k)]; count++; }
                    if (j > 1) { sum += u_new[U_IX(i, j-1, k)]; count++; }
                    if (j < N) { sum += u_new[U_IX(i, j+1, k)]; count++; }
                    if (k > 1) { sum += u_new[U_IX(i, j, k-1)]; count++; }
                    if (k < N) { sum += u_new[U_IX(i, j, k+1)]; count++; }

                    if (count > 0) {
                        // Use u_rhs (advected velocities) as right-hand side
                        // Equation: (I - ν∆t∇²)u_new = u_rhs
                        // Discretized: (1 + a*count)*u_new - a*sum = u_rhs
                        // Rearranged: u_new = (u_rhs + a*sum) / (1 + a*count)
                        u_new[idx] = (u_rhs[idx] + a * sum) / (1.0f + a * count);
                    }
                }
            }
        }
        set_bnd(1, u_new, N+1, N+2, N+2);
    }
    u = std::move(u_new);

    // Solve for v (similar logic)
    for (int iter = 0; iter < 20; iter++) {
        for (int k = 1; k <= N; k++) {
            for (int j = 0; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
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
        set_bnd(2, v_new, N+2, N+1, N+2);
    }
    v = std::move(v_new);

    // Solve for w (similar logic)
    for (int iter = 0; iter < 20; iter++) {
        for (int k = 0; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
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
        set_bnd(3, w_new, N+2, N+2, N+1);
    }
    w = std::move(w_new);
}

    void diffuse_density(float dt) {
        if (diff <= 0.0f) return;

        float a = dt * diff * N * N;

        // Save advected density as right-hand side
        std::vector<float> dens_rhs = dens;
        std::vector<float> dens_new = dens;

        for (int iter = 0; iter < 20; iter++) {
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        int idx = P_IX(i, j, k);
                        // Use dens_new for neighbors (Gauss-Seidel)
                        float sum = dens_new[P_IX(i-1, j, k)] + dens_new[P_IX(i+1, j, k)] +
                                    dens_new[P_IX(i, j-1, k)] + dens_new[P_IX(i, j+1, k)] +
                                    dens_new[P_IX(i, j, k-1)] + dens_new[P_IX(i, j, k+1)];

                        // Use dens_rhs (advected density) as right-hand side
                        dens_new[idx] = (dens_rhs[idx] + a * sum) / (1 + 6 * a);
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

        return (u_right - u_left) + (v_top - v_bottom) + (w_front - w_back);
    }

    void project() {

        int size = N + 2;

        // 1. Compute divergence of velocity field w3
        std::vector<float> div((N+2)*(N+2)*(N+2), 0.0f);
        float max_div_before = 0.0f;

        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    float d = computeDivergence(i, j, k);
                    div[P_IX(i, j, k)] = d;
                    if (fabs(d) > max_div_before) max_div_before = fabs(d);
                }
            }
        }
        set_bnd(0, div, size, size, size);

        std::cout << "[PROJECT] Max divergence before: " << max_div_before << std::endl;

        // 2. Solve Poisson equation: ∇²p = ∇·u
        // Actually, from paper: ∇²q = ∇·w3
        // Where q is pressure-like scalar field

        // Clear pressure field each timestep
        std::fill(p.begin(), p.end(), 0.0f);

        // SOR parameters
        const int max_iterations = 20;  // Increased for better convergence
        const float tolerance = 1e-4f;
        const float omega = 1.85f;

        // CORRECT SCALE: In Poisson equation ∇²p = ∇·u, there's no scale!
        // The RHS is just divergence, not scaled divergence
        // So scale = 1.0f, NOT rho*dx²/dt

        for (int iter = 0; iter < max_iterations; iter++) {
            float max_change = 0.0f;

            // Red cells
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    int start_i = ((j + k) % 2 == 0) ? 1 : 2;
                    for (int i = start_i; i <= N; i += 2) {
                        int idx = P_IX(i, j, k);

                        float sum_neighbors = p[P_IX(i-1, j, k)] + p[P_IX(i+1, j, k)] +
                                             p[P_IX(i, j-1, k)] + p[P_IX(i, j+1, k)] +
                                             p[P_IX(i, j, k-1)] + p[P_IX(i, j, k+1)];

                        // RHS is just divergence, no scaling!
                        float rhs = div[idx];

                        // Gauss-Seidel with SOR: p_new = (sum_neighbors - rhs)/6
                        // But for Poisson: sum_neighbors - 6*p = rhs
                        // So: p = (sum_neighbors - rhs)/6
                        float p_new = (1.0f - omega) * p[idx] +
                                     omega * (sum_neighbors - rhs) / 6.0f;

                        float change = fabs(p_new - p[idx]);
                        if (change > max_change) max_change = change;

                        p[idx] = p_new;
                    }
                }
            }
            set_bnd(0, p, size, size, size);

            // Black cells
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    int start_i = ((j + k) % 2 == 0) ? 2 : 1;
                    for (int i = start_i; i <= N; i += 2) {
                        int idx = P_IX(i, j, k);

                        float sum_neighbors = p[P_IX(i-1, j, k)] + p[P_IX(i+1, j, k)] +
                                             p[P_IX(i, j-1, k)] + p[P_IX(i, j+1, k)] +
                                             p[P_IX(i, j, k-1)] + p[P_IX(i, j, k+1)];

                        float rhs = div[idx];
                        float p_new = (1.0f - omega) * p[idx] +
                                     omega * (sum_neighbors - rhs) / 6.0f;

                        float change = fabs(p_new - p[idx]);
                        if (change > max_change) max_change = change;

                        p[idx] = p_new;
                    }
                }
            }
            set_bnd(0, p, size, size, size);

            if (max_change < tolerance) {
                std::cout << "[PROJECT] Converged after " << iter+1 << " iterations" << std::endl;
                break;
            }
    }

    // 3. Apply pressure gradient: w4 = w3 - ∇p
    applyPressureGradient();

    // Check divergence after projection
    float max_div_after = 0.0f;
    for (int k = 1; k <= N; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = 1; i <= N; i++) {
                float d = computeDivergence(i, j, k);
                if (fabs(d) > max_div_after) max_div_after = fabs(d);
            }
        }
    }

    std::cout << "[PROJECT] Max divergence after: " << max_div_after << std::endl;
    std::cout << "[PROJECT] Divergence reduced by factor: "
              << (max_div_before > 0 ? max_div_after / max_div_before : 0) << std::endl;
    }

    void applyPressureGradient() {
        // Apply gradient: w4 = w3 - ∇p
        // Note: No dt scaling here! The pressure p already incorporates any scaling
        // from solving the Poisson equation

        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 0; i <= N; i++) {
                    float p_right = p[P_IX(i+1, j, k)];
                    float p_left = p[P_IX(i, j, k)];
                    float pressure_grad = p_right - p_left;

                    u[U_IX(i, j, k)] -= pressure_grad;  // No dt scaling!
                }
            }
        }

        for (int k = 1; k <= N; k++) {
            for (int j = 0; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    float p_top = p[P_IX(i, j+1, k)];
                    float p_bottom = p[P_IX(i, j, k)];
                    float pressure_grad = p_top - p_bottom;

                    v[V_IX(i, j, k)] -= pressure_grad;  // No dt scaling!
                }
            }
        }

        for (int k = 0; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    float p_front = p[P_IX(i, j, k+1)];
                    float p_back = p[P_IX(i, j, k)];
                    float pressure_grad = p_front - p_back;

                    w[W_IX(i, j, k)] -= pressure_grad;  // No dt scaling!
                }
            }
        }

        // Apply boundary conditions
        set_bnd(1, u, N+1, N+2, N+2);
        set_bnd(2, v, N+2, N+1, N+2);
        set_bnd(3, w, N+2, N+2, N+1);
    }

    void dissipate_density(float dt, float alpha = 0.1f) {
        // Dissipation: d_new = d_old / (1 + ∆tα)
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
        // Add density at bottom center
        int centerX = N / 2;
        int centerY = 2;
        int centerZ = N / 2;

        for (int k = centerZ-1; k <= centerZ+1; k++) {
            for (int i = centerX-1; i <= centerX+1; i++) {
                add_density(i, centerY, k, 150.0f);
            }
        }

        // Add upward velocity
        for (int k = centerZ-1; k <= centerZ+1; k++) {
            for (int i = centerX-1; i <= centerX+1; i++) {
                add_velocity(i, centerY, k, 0.0f, 1.5f, 0.0f);
            }
        }
    }

    void vel_step() {
        add_forces();
        advect_velocity(dt);
        diffuse_velocity(dt);
        project();
    }

    void dens_step() {
        advect_density(dt);
        diffuse_density(dt);
    }


};

#endif