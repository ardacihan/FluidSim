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

private:

    bool is_solid(int x, int y, int z) {
        if (x==0 || y==0 || z==0 || x==N || y==N || z==N) {
            return true;
        }
        return false;
    }

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
        // Copy current velocities to old arrays
        std::copy(u.begin(), u.end(), u_old.begin());
        std::copy(v.begin(), v.end(), v_old.begin());
        std::copy(w.begin(), w.end(), w_old.begin());

        // Advect u
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 0; i <= N; i++) {
                    float x = (float)i + 0.5f;
                    float y = (float)j + 0.5f;
                    float z = (float)k + 0.5f;

                    float u_vel = u_old[U_IX(i, j, k)];
                    float v_vel = 0.25f * (
                        v_old[V_IX(i, j-1, k)] + v_old[V_IX(i, j, k)] +
                        v_old[V_IX(i+1, j-1, k)] + v_old[V_IX(i+1, j, k)]
                    );
                    float w_vel = 0.25f * (
                        w_old[W_IX(i, j, k-1)] + w_old[W_IX(i, j, k)] +
                        w_old[W_IX(i+1, j, k-1)] + w_old[W_IX(i+1, j, k)]
                    );

                    float srcX = x - dt * u_vel;
                    float srcY = y - dt * v_vel;
                    float srcZ = z - dt * w_vel;

                    u[U_IX(i, j, k)] = interpolate_u(srcX, srcY, srcZ);
                }
            }
        }
         set_bnd(1, u, N+1, N+2, N+2);

        //v
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

                    float srcX = x - dt * u_vel;
                    float srcY = y - dt * v_vel;
                    float srcZ = z - dt * w_vel;

                    v[V_IX(i, j, k)] = interpolate_v(srcX, srcY, srcZ);
                }
            }
        }
        // boundary  v
        set_bnd(2, v, N+2, N+1, N+2);

        // w
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

                    float srcX = x - dt * u_vel;
                    float srcY = y - dt * v_vel;
                    float srcZ = z - dt * w_vel;

                    w[W_IX(i, j, k)] = interpolate_w(srcX, srcY, srcZ);
                }
            }
        }
        // boundary w
        set_bnd(3, w, N+2, N+2, N+1);
    }

    void advect_density(float dt) {

    }

    void vel_step() {
        add_velocity(0,0,0,0,0,0);           // Add gravity, user forces, etc.
        diffuse_velocity(dt);     // Apply viscosity (if visc > 0)
        project();               // Make divergence-free (pressure solve)
        advect_velocity(dt);     // Move velocity with itself
        project();               // Make divergence-free again
    }

    void dens_step() {
        add_density(0,0,0,0);  // Add new density
        diffuse_density(dt);      // Apply diffusion (if diff > 0)
        advect_density(dt);       // Move density with velocity
    }

    void add_source(int u,int v, int w, int size_x, int size_y, int size_z) {
    }

    void diffuse_velocity(float dt) {
    }
    void diffuse_density(float dt) {
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
        float rho = 1.0f;  // Density (water = 1000 kg/m³, but we use 1.0)
        float dx = 1.0f;   // Grid spacing

        // Step 1: Compute divergence at each fluid cell
        std::vector<float> div((N+2)*(N+2)*(N+2), 0.0f);
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    div[P_IX(i, j, k)] = computeDivergence(i, j, k);
                }
            }
        }
        set_bnd(0, div, size, size, size);

        // Step 2: Solve the linear system: A * p = b
        // Where A is: 6p - sum_neighbors = (ρ * dx² / Δt) * (-div)
        // Rearranged: p = (sum_neighbors - (ρ * dx² / Δt) * div) / 6

        float scale = rho * dx * dx / dt;
        std::vector<float> p_new = p;  // Copy current pressure

        for (int iter = 0; iter < 50; iter++) {
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        int idx = P_IX(i, j, k);

                        // Get neighbor pressures
                        float sum_neighbors = 0.0f;
                        int neighbor_count = 0;

                        // Check each neighbor (fluid cells only)
                        if (i > 1) {
                            sum_neighbors += p[P_IX(i-1, j, k)];
                            neighbor_count++;
                        }
                        if (i < N) {
                            sum_neighbors += p[P_IX(i+1, j, k)];
                            neighbor_count++;
                        }
                        if (j > 1) {
                            sum_neighbors += p[P_IX(i, j-1, k)];
                            neighbor_count++;
                        }
                        if (j < N) {
                            sum_neighbors += p[P_IX(i, j+1, k)];
                            neighbor_count++;
                        }
                        if (k > 1) {
                            sum_neighbors += p[P_IX(i, j, k-1)];
                            neighbor_count++;
                        }
                        if (k < N) {
                            sum_neighbors += p[P_IX(i, j, k+1)];
                            neighbor_count++;
                        }

                        // Apply equation from paper
                        if (neighbor_count == 6) {
                            // Interior cell: all 6 neighbors are fluid
                            p_new[idx] = (sum_neighbors - scale * div[idx]) / 6.0f;
                        } else if (neighbor_count > 0) {
                            // Boundary cell: some neighbors are solid
                            p_new[idx] = (sum_neighbors - scale * div[idx]) / (float)neighbor_count;
                        }
                    }
                }
            }

            // Swap p and p_new
            std::swap(p, p_new);
            set_bnd(0, p, size, size, size);
        }

        // Step 3: Apply pressure gradient to velocities
        applyPressureGradient();
    }


    void applyPressureGradient() {
        float pressure_scale = dt;  // dt/ρ with ρ=1.0

        // Update u velocities
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 0; i <= N; i++) {
                    // Get velocity at this u-face BEFORE pressure
                    auto vel = getVelocityAtUFace(i, j, k);
                    float u_before = vel[0];

                    // Pressure gradient: (p_right - p_left)
                    float p_right = p[P_IX(i+1, j, k)];
                    float p_left = p[P_IX(i, j, k)];
                    float pressure_grad = p_right - p_left;

                    // Update: u_new = u_old - (Δt/ρ) * ∇p
                    u[U_IX(i, j, k)] = u_before - pressure_scale * pressure_grad;
                }
            }
        }

        // Update v velocities (y-faces)
        for (int k = 1; k <= N; k++) {
            for (int j = 0; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    auto vel = getVelocityAtVFace(i, j, k);
                    float v_before = vel[1];

                    float p_top = p[P_IX(i, j+1, k)];
                    float p_bottom = p[P_IX(i, j, k)];
                    float pressure_grad = p_top - p_bottom;

                    v[V_IX(i, j, k)] = v_before - pressure_scale * pressure_grad;
                }
            }
        }

        // Update w velocities
        for (int k = 0; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    auto vel = getVelocityAtWFace(i, j, k);
                    float w_before = vel[2];

                    float p_front = p[P_IX(i, j, k+1)];
                    float p_back = p[P_IX(i, j, k)];
                    float pressure_grad = p_front - p_back;

                    w[W_IX(i, j, k)] = w_before - pressure_scale * pressure_grad;
                }
            }
        }

        // Apply boundary conditions
        set_bnd(1, u, N+1, N+2, N+2);
        set_bnd(2, v, N+2, N+1, N+2);
        set_bnd(3, w, N+2, N+2, N+1);
    }
};

#endif