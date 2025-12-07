#ifndef GRID3D_H
#define GRID3D_H

#include <vector>
#include <algorithm>
#include <cmath>

class Grid3D {
public:
    int N; // Grid size (N x N x N) - number of pressure cells
    float dt;
    float diff;
    float visc;

    // Pressure at cell center
    std::vector<float> p;

    // Velocities (x,y,z) (MAC/ staggered grid)
    std::vector<float> u;
    std::vector<float> v;
    std::vector<float> w;

    std::vector<float> u_old, v_old, w_old;

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

private:

    bool is_solid(int x, int y, int z) {
        if (x==0 || y==0 || z==0 || x==N || y==N || z==N) {
            return true;
        }
        return false;

    }

    void set_bnd(int b, std::vector<float>& x, int size_x, int size_y, int size_z) {
    }


    // Interpolate density (cell-centered quantity)
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

    float interpolate_u(float x, float y, float z) const {
        x = std::max(0.5f, std::min((float)N + 0.5f, x));
        y = std::max(0.5f, std::min((float)N + 1.5f, y));
        z = std::max(0.5f, std::min((float)N + 1.5f, z));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);
        int k = (int)floor(z - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;
        float u = (z - 0.5f) - k;

        i = std::max(0, std::min(N, i));
        j = std::max(0, std::min(N+1, j));
        k = std::max(0, std::min(N+1, k));

        return (1-s)*(1-t)*(1-u)*u_old[U_IX(i, j, k)] +
               s*(1-t)*(1-u)*u_old[U_IX(i+1, j, k)] +
               (1-s)*t*(1-u)*u_old[U_IX(i, j+1, k)] +
               s*t*(1-u)*u_old[U_IX(i+1, j+1, k)] +
               (1-s)*(1-t)*u*u_old[U_IX(i, j, k+1)] +
               s*(1-t)*u*u_old[U_IX(i+1, j, k+1)] +
               (1-s)*t*u*u_old[U_IX(i, j+1, k+1)] +
               s*t*u*u_old[U_IX(i+1, j+1, k+1)];
    }

    // Interpolate v (y-face quantity)
    float interpolate_v(float x, float y, float z) const {
        x = std::max(0.5f, std::min((float)N + 1.5f, x));
        y = std::max(0.5f, std::min((float)N + 0.5f, y));
        z = std::max(0.5f, std::min((float)N + 1.5f, z));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);
        int k = (int)floor(z - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;
        float u = (z - 0.5f) - k;

        i = std::max(0, std::min(N+1, i));
        j = std::max(0, std::min(N, j));
        k = std::max(0, std::min(N+1, k));

        return (1-s)*(1-t)*(1-u)*v_old[V_IX(i, j, k)] +
               s*(1-t)*(1-u)*v_old[V_IX(i+1, j, k)] +
               (1-s)*t*(1-u)*v_old[V_IX(i, j+1, k)] +
               s*t*(1-u)*v_old[V_IX(i+1, j+1, k)] +
               (1-s)*(1-t)*u*v_old[V_IX(i, j, k+1)] +
               s*(1-t)*u*v_old[V_IX(i+1, j, k+1)] +
               (1-s)*t*u*v_old[V_IX(i, j+1, k+1)] +
               s*t*u*v_old[V_IX(i+1, j+1, k+1)];
    }

    // Interpolate w (z-face quantity)
    float interpolate_w(float x, float y, float z) const {
        x = std::max(0.5f, std::min((float)N + 1.5f, x));
        y = std::max(0.5f, std::min((float)N + 1.5f, y));
        z = std::max(0.5f, std::min((float)N + 0.5f, z));

        int i = (int)floor(x - 0.5f);
        int j = (int)floor(y - 0.5f);
        int k = (int)floor(z - 0.5f);

        float s = (x - 0.5f) - i;
        float t = (y - 0.5f) - j;
        float u = (z - 0.5f) - k;

        i = std::max(0, std::min(N+1, i));
        j = std::max(0, std::min(N+1, j));
        k = std::max(0, std::min(N, k));

        return (1-s)*(1-t)*(1-u)*w_old[W_IX(i, j, k)] +
               s*(1-t)*(1-u)*w_old[W_IX(i+1, j, k)] +
               (1-s)*t*(1-u)*w_old[W_IX(i, j+1, k)] +
               s*t*(1-u)*w_old[W_IX(i+1, j+1, k)] +
               (1-s)*(1-t)*u*w_old[W_IX(i, j, k+1)] +
               s*(1-t)*u*w_old[W_IX(i+1, j, k+1)] +
               (1-s)*t*u*w_old[W_IX(i, j+1, k+1)] +
               s*t*u*w_old[W_IX(i+1, j+1, k+1)];
    }


    void project() {
    }

    void advect(std::vector<float>& d, std::vector<float>& d0,
                const std::vector<float>& u, const std::vector<float>& v, const std::vector<float>& w,
                int size_x, int size_y, int size_z, float dt) {
    }

    void vel_step() {
    }

    void dens_step() {

    }

    void add_source(std::vector<float>& x, const std::vector<float>& s, float dt,
                    int size_x, int size_y, int size_z) {
    }

    void diffuse(int b, std::vector<float>& x, std::vector<float>& x0, float diff, float dt,
                 int size_x, int size_y, int size_z) {
    }
};

#endif