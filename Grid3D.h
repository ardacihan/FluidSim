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
    int N; // Grid size (N x N x N)
    float dt;
    float diff;
    float visc;
    float h = 1.0f;
    float inv_h = 1.0f;

    std::vector<float> p;
    std::vector<float> u, v, w;
    std::vector<float> u_old, v_old, w_old;
    std::vector<float> dens, dens_old;

    Grid3D(int size, float diffusion, float viscosity, float timestep)
        : N(size), diff(diffusion), visc(viscosity), dt(timestep)
    {
        int cell_count = (N + 2) * (N + 2) * (N + 2);
        p.resize(cell_count, 0.0f);
        dens.resize(cell_count, 0.0f);
        dens_old.resize(cell_count, 0.0f);

        int u_count = (N + 1) * (N + 2) * (N + 2);
        u.resize(u_count, 0.0f);
        u_old.resize(u_count, 0.0f);

        int v_count = (N + 2) * (N + 1) * (N + 2);
        v.resize(v_count, 0.0f);
        v_old.resize(v_count, 0.0f);

        int w_count = (N + 2) * (N + 2) * (N + 1);
        w.resize(w_count, 0.0f);
        w_old.resize(w_count, 0.0f);

        h = 1.0f / N;
        inv_h = N;
    }

    inline int P_IX(int i, int j, int k) const {
        return i + (N + 2) * j + (N + 2) * (N + 2) * k;
    }

    inline int U_IX(int i, int j, int k) const {
        return i + (N + 1) * j + (N + 1) * (N + 2) * k;
    }

    inline int V_IX(int i, int j, int k) const {
        return i + (N + 2) * j + (N + 2) * (N + 1) * k;
    }

    inline int W_IX(int i, int j, int k) const {
        return i + (N + 2) * j + (N + 2) * (N + 2) * k;
    }

    void add_density(int i, int j, int k, float amount) {
        dens[P_IX(i, j, k)] += amount;
    }

    void add_velocity(int i, int j, int k, float amountU, float amountV, float amountW) {
        u[U_IX(i, j, k)] += amountU * 0.5f;
        u[U_IX(i+1, j, k)] += amountU * 0.5f;
        v[V_IX(i, j, k)] += amountV * 0.5f;
        v[V_IX(i, j+1, k)] += amountV * 0.5f;
        w[W_IX(i, j, k)] += amountW * 0.5f;
        w[W_IX(i, j, k+1)] += amountW * 0.5f;
    }

    void step() {
        vel_step();
        dens_step();
    }

    std::vector<float> getVelocityAtCellCenter(int i, int j, int k) const {
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

    void initRandomStaggeredVelocities(float magnitude = 0.5f) {
        clearAllVelocities();

        for (int k = 0; k <= N+1; k++) {
            for (int j = 0; j <= N+1; j++) {
                for (int i = 0; i <= N; i++) {
                    u[U_IX(i, j, k)] = magnitude * ((rand() / (float)RAND_MAX) * 2.0f - 1.0f);
                }
            }
        }

        for (int k = 0; k <= N+1; k++) {
            for (int j = 0; j <= N; j++) {
                for (int i = 0; i <= N+1; i++) {
                    v[V_IX(i, j, k)] = magnitude * ((rand() / (float)RAND_MAX) * 2.0f - 1.0f);
                }
            }
        }

        for (int k = 0; k <= N; k++) {
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
        #pragma omp parallel for collapse(2)
        for (int k = 1; k < size_z-1; k++) {
            for (int j = 1; j < size_y-1; j++) {
                if (b == 1) {
                    x[0 + j*size_x + k*size_x*size_y] = 0.0f;
                    x[(size_x-1) + j*size_x + k*size_x*size_y] = 0.0f;
                } else {
                    x[0 + j*size_x + k*size_x*size_y] = x[1 + j*size_x + k*size_x*size_y];
                    x[(size_x-1) + j*size_x + k*size_x*size_y] = x[(size_x-2) + j*size_x + k*size_x*size_y];
                }
            }
        }

        #pragma omp parallel for collapse(2)
        for (int k = 1; k < size_z-1; k++) {
            for (int i = 1; i < size_x-1; i++) {
                if (b == 2) {
                    x[i + 0*size_x + k*size_x*size_y] = 0.0f;
                    x[i + (size_y-1)*size_x + k*size_x*size_y] = 0.0f;
                } else {
                    x[i + 0*size_x + k*size_x*size_y] = x[i + 1*size_x + k*size_x*size_y];
                    x[i + (size_y-1)*size_x + k*size_x*size_y] = x[i + (size_y-2)*size_x + k*size_x*size_y];
                }
            }
        }

        #pragma omp parallel for collapse(2)
        for (int j = 1; j < size_y-1; j++) {
            for (int i = 1; i < size_x-1; i++) {
                if (b == 3) {
                    x[i + j*size_x + 0*size_x*size_y] = 0.0f;
                    x[i + j*size_x + (size_z-1)*size_x*size_y] = 0.0f;
                } else {
                    x[i + j*size_x + 0*size_x*size_y] = x[i + j*size_x + 1*size_x*size_y];
                    x[i + j*size_x + (size_z-1)*size_x*size_y] = x[i + j*size_x + (size_z-2)*size_x*size_y];
                }
            }
        }
    }

    void advect_velocity(float dt) {
        std::copy(u.begin(), u.end(), u_old.begin());
        std::copy(v.begin(), v.end(), v_old.begin());
        std::copy(w.begin(), w.end(), w_old.begin());

        // Simplified backward Euler advection (faster than RK2)
        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 0; i <= N; i++) {
                    float x = i + 0.5f;
                    float y = j + 0.5f;
                    float z = k + 0.5f;

                    float u_vel = u_old[U_IX(i, j, k)];
                    float v_vel = 0.25f * (
                        v_old[V_IX(i, j-1, k)] + v_old[V_IX(i, j, k)] +
                        v_old[V_IX(i+1, j-1, k)] + v_old[V_IX(i+1, j, k)]
                    );
                    float w_vel = 0.25f * (
                        w_old[W_IX(i, j, k-1)] + w_old[W_IX(i, j, k)] +
                        w_old[W_IX(i+1, j, k-1)] + w_old[W_IX(i+1, j, k)]
                    );

                    float src_x = x - dt * u_vel * inv_h;
                    float src_y = y - dt * v_vel * inv_h;
                    float src_z = z - dt * w_vel * inv_h;

                    u[U_IX(i, j, k)] = interpolate_u_old(src_x, src_y, src_z);
                }
            }
        }

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

                    float src_x = x - dt * u_vel * inv_h;
                    float src_y = y - dt * v_vel * inv_h;
                    float src_z = z - dt * w_vel * inv_h;

                    v[V_IX(i, j, k)] = interpolate_v_old(src_x, src_y, src_z);
                }
            }
        }

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

                    float src_x = x - dt * u_vel * inv_h;
                    float src_y = y - dt * v_vel * inv_h;
                    float src_z = z - dt * w_vel * inv_h;

                    w[W_IX(i, j, k)] = interpolate_w_old(src_x, src_y, src_z);
                }
            }
        }

        set_bnd(1, u, N+1, N+2, N+2);
        set_bnd(2, v, N+2, N+1, N+2);
        set_bnd(3, w, N+2, N+2, N+1);
    }

    void advect_density(float dt) {
        std::copy(dens.begin(), dens.end(), dens_old.begin());

        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    float x = (float)i + 0.5f;
                    float y = (float)j + 0.5f;
                    float z = (float)k + 0.5f;

                    auto vel = getVelocityAtCellCenter(i, j, k);
                    float u_vel = vel[0];
                    float v_vel = vel[1];
                    float w_vel = vel[2];

                    float srcX = x - dt * u_vel * inv_h;
                    float srcY = y - dt * v_vel * inv_h;
                    float srcZ = z - dt * w_vel * inv_h;

                    srcX = std::max(0.5f, std::min((float)N + 0.5f, srcX));
                    srcY = std::max(0.5f, std::min((float)N + 0.5f, srcY));
                    srcZ = std::max(0.5f, std::min((float)N + 0.5f, srcZ));

                    dens[P_IX(i, j, k)] = interpolate_density(srcX, srcY, srcZ);
                }
            }
        }
        set_bnd(0, dens, N+2, N+2, N+2);
    }

    void diffuse_velocity(float dt) {
        if (visc <= 0.0f) return;

        float a = dt * visc * inv_h * inv_h;
        float factor = 1.0f / (1.0f + 6.0f * a);

        std::vector<float> u_temp = u;
        std::vector<float> v_temp = v;
        std::vector<float> w_temp = w;

        // Reduced iterations with Jacobi (better for parallel)
        for (int iter = 0; iter < 5; iter++) {
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 0; i <= N; i++) {
                        int idx = U_IX(i, j, k);
                        float sum = 0.0f;
                        int count = 0;

                        if (i > 0) { sum += u_temp[U_IX(i-1, j, k)]; count++; }
                        if (i < N) { sum += u_temp[U_IX(i+1, j, k)]; count++; }
                        if (j > 1) { sum += u_temp[U_IX(i, j-1, k)]; count++; }
                        if (j < N) { sum += u_temp[U_IX(i, j+1, k)]; count++; }
                        if (k > 1) { sum += u_temp[U_IX(i, j, k-1)]; count++; }
                        if (k < N) { sum += u_temp[U_IX(i, j, k+1)]; count++; }

                        if (count > 0) {
                            u[idx] = (u_temp[idx] + a * sum) / (1.0f + a * count);
                        }
                    }
                }
            }
            std::swap(u, u_temp);

            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 0; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        int idx = V_IX(i, j, k);
                        float sum = 0.0f;
                        int count = 0;

                        if (i > 1) { sum += v_temp[V_IX(i-1, j, k)]; count++; }
                        if (i < N) { sum += v_temp[V_IX(i+1, j, k)]; count++; }
                        if (j > 0) { sum += v_temp[V_IX(i, j-1, k)]; count++; }
                        if (j < N) { sum += v_temp[V_IX(i, j+1, k)]; count++; }
                        if (k > 1) { sum += v_temp[V_IX(i, j, k-1)]; count++; }
                        if (k < N) { sum += v_temp[V_IX(i, j, k+1)]; count++; }

                        if (count > 0) {
                            v[idx] = (v_temp[idx] + a * sum) / (1.0f + a * count);
                        }
                    }
                }
            }
            std::swap(v, v_temp);

            #pragma omp parallel for collapse(3)
            for (int k = 0; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        int idx = W_IX(i, j, k);
                        float sum = 0.0f;
                        int count = 0;

                        if (i > 1) { sum += w_temp[W_IX(i-1, j, k)]; count++; }
                        if (i < N) { sum += w_temp[W_IX(i+1, j, k)]; count++; }
                        if (j > 1) { sum += w_temp[W_IX(i, j-1, k)]; count++; }
                        if (j < N) { sum += w_temp[W_IX(i, j+1, k)]; count++; }
                        if (k > 0) { sum += w_temp[W_IX(i, j, k-1)]; count++; }
                        if (k < N) { sum += w_temp[W_IX(i, j, k+1)]; count++; }

                        if (count > 0) {
                            w[idx] = (w_temp[idx] + a * sum) / (1.0f + a * count);
                        }
                    }
                }
            }
            std::swap(w, w_temp);
        }

        u = u_temp;
        v = v_temp;
        w = w_temp;

        set_bnd(1, u, N+1, N+2, N+2);
        set_bnd(2, v, N+2, N+1, N+2);
        set_bnd(3, w, N+2, N+2, N+1);
    }

    void diffuse_density(float dt) {
        if (diff <= 0.0f) return;

        float a = dt * diff * inv_h * inv_h;
        float factor = 1.0f / (1.0f + 6.0f * a);

        std::vector<float> dens_temp = dens;

        for (int iter = 0; iter < 10; iter++) {
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        int idx = P_IX(i, j, k);
                        float sum = dens_temp[P_IX(i-1, j, k)] + dens_temp[P_IX(i+1, j, k)] +
                                    dens_temp[P_IX(i, j-1, k)] + dens_temp[P_IX(i, j+1, k)] +
                                    dens_temp[P_IX(i, j, k-1)] + dens_temp[P_IX(i, j, k+1)];

                        dens[idx] = (dens_temp[idx] + a * sum) * factor;
                    }
                }
            }
            std::swap(dens, dens_temp);
        }
        dens = dens_temp;
        set_bnd(0, dens, N+2, N+2, N+2);
    }

    void project(float dt) {
        std::vector<float> div((N+2)*(N+2)*(N+2), 0.0f);
        std::vector<float> p((N+2)*(N+2)*(N+2), 0.0f);
        std::vector<float> p_temp = p;

        // Compute divergence
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

        float h_squared = h * h;
        float inv_6 = 1.0f / 6.0f;

        // Solve Poisson equation with Jacobi (15 iterations for balance)
        for (int iter = 0; iter < 15; iter++) {
            #pragma omp parallel for collapse(3)
            for (int k = 1; k <= N; k++) {
                for (int j = 1; j <= N; j++) {
                    for (int i = 1; i <= N; i++) {
                        int idx = P_IX(i, j, k);
                        float p_sum = p_temp[P_IX(i-1, j, k)] + p_temp[P_IX(i+1, j, k)] +
                                      p_temp[P_IX(i, j-1, k)] + p_temp[P_IX(i, j+1, k)] +
                                      p_temp[P_IX(i, j, k-1)] + p_temp[P_IX(i, j, k+1)];

                        p[idx] = (p_sum - h_squared * div[idx]) * inv_6;
                    }
                }
            }
            std::swap(p, p_temp);

            // Update boundaries every 5 iterations
            if (iter % 5 == 0) {
                set_bnd(0, p_temp, N+2, N+2, N+2);
            }
        }
        p = p_temp;
        set_bnd(0, p, N+2, N+2, N+2);

        // Apply pressure gradient
        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N-1; i++) {
                    u[U_IX(i, j, k)] -= dt * (p[P_IX(i+1, j, k)] - p[P_IX(i, j, k)]) * inv_h;
                }
            }
        }

        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N; k++) {
            for (int j = 1; j <= N-1; j++) {
                for (int i = 1; i <= N; i++) {
                    v[V_IX(i, j, k)] -= dt * (p[P_IX(i, j+1, k)] - p[P_IX(i, j, k)]) * inv_h;
                }
            }
        }

        #pragma omp parallel for collapse(3)
        for (int k = 1; k <= N-1; k++) {
            for (int j = 1; j <= N; j++) {
                for (int i = 1; i <= N; i++) {
                    w[W_IX(i, j, k)] -= dt * (p[P_IX(i, j, k+1)] - p[P_IX(i, j, k)]) * inv_h;
                }
            }
        }

        set_bnd(1, u, N+1, N+2, N+2);
        set_bnd(2, v, N+2, N+1, N+2);
        set_bnd(3, w, N+2, N+2, N+1);

        this->p = p;
    }

    void dissipate_density(float dt, float alpha = 0.1f) {
        float factor = 1.0f / (1.0f + dt * alpha);

        #pragma omp parallel for collapse(3)
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

        for (int k = centerZ-2; k <= centerZ+2; k++) {
            for (int i = centerX-2; i <= centerX+2; i++) {
                float dx = (i - centerX) * 0.3f;
                float dz = (k - centerZ) * 0.3f;
                add_velocity(i, centerY, k, -dz, 1.5f, dx);
                add_density(i, centerY, k, 150.0f);
            }
        }
    }

    void vel_step() {
        add_forces();
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