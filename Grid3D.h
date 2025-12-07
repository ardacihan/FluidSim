#ifndef GRID3D_H
#define GRID3D_H

#include <vector>
#include <algorithm>
#include <cmath>

class Grid3D {
public:
    int N; // Grid size (N x N x N)
    float dt;
    float diff;
    float visc;

    // Flattened arrays for 3D data
    // size = (N+2)*(N+2)*(N+2)
    std::vector<float> u, v, w;     // Velocity (x, y, z)
    std::vector<float> u0, v0, w0;  // Previous velocity
    std::vector<float> dens;        // Density
    std::vector<float> dens0;       // Previous density

    Grid3D(int size, float diffusion, float viscosity, float timestep)
        : N(size), diff(diffusion), visc(viscosity), dt(timestep)
    {
        int total_size = (N + 2) * (N + 2) * (N + 2);
        u.resize(total_size, 0.0f); v.resize(total_size, 0.0f); w.resize(total_size, 0.0f);
        u0.resize(total_size, 0.0f); v0.resize(total_size, 0.0f); w0.resize(total_size, 0.0f);
        dens.resize(total_size, 0.0f); dens0.resize(total_size, 0.0f);
    }

    // Helper to turn (x,y,z) into 1D index
    // IX(x, y, z) = x + (N+2)*y + (N+2)*(N+2)*z
    inline int IX(int x, int y, int z) const {
        return x + (N + 2) * y + (N + 2) * (N + 2) * z;
    }

    void add_density(int x, int y, int z, float amount) {
        dens[IX(x, y, z)] += amount;
    }

    void add_velocity(int x, int y, int z, float amountX, float amountY, float amountZ) {
        int index = IX(x, y, z);
        u[index] += amountX;
        v[index] += amountY;
        w[index] += amountZ;
    }

    void step() {
        vel_step(u, v, w, u0, v0, w0, visc, dt);
        dens_step(dens, dens0, u, v, w, diff, dt);
    }

private:
    void set_bnd(int b, std::vector<float>& x) {
        //walls reflect
        for(int j=1; j<=N; j++) {
            for(int i=1; i<=N; i++) {
                x[IX(i, j, 0  )] = b==3 ? -x[IX(i, j, 1  )] : x[IX(i, j, 1  )];
                x[IX(i, j, N+1)] = b==3 ? -x[IX(i, j, N  )] : x[IX(i, j, N  )];
            }
        }
        for(int k=1; k<=N; k++) {
            for(int i=1; i<=N; i++) {
                x[IX(i, 0, k  )] = b==2 ? -x[IX(i, 1, k  )] : x[IX(i, 1, k  )];
                x[IX(i, N+1, k)] = b==2 ? -x[IX(i, N, k  )] : x[IX(i, N, k  )];
            }
        }
        for(int k=1; k<=N; k++) {
            for(int j=1; j<=N; j++) {
                x[IX(0, j, k  )] = b==1 ? -x[IX(1, j, k  )] : x[IX(1, j, k  )];
                x[IX(N+1, j, k)] = b==1 ? -x[IX(N, j, k  )] : x[IX(N, j, k  )];
            }
        }
        x[IX(0,0,0)] = 0.33f * (x[IX(1,0,0)] + x[IX(0,1,0)] + x[IX(0,0,1)]);
    }

    void lin_solve(int b, std::vector<float>& x, std::vector<float>& x0, float a, float c) {
    }

    void project(std::vector<float>& velocX, std::vector<float>& velocY, std::vector<float>& velocZ, std::vector<float>& p, std::vector<float>& div) {
    }

    void advect(int b, std::vector<float>& d, std::vector<float>& d0, std::vector<float>& velocX, std::vector<float>& velocY, std::vector<float>& velocZ, float dt) {
    }

    void vel_step(std::vector<float>& x, std::vector<float>& y, std::vector<float>& z, std::vector<float>& x0, std::vector<float>& y0, std::vector<float>& z0, float visc, float dt) {
    }

    void dens_step(std::vector<float>& x, std::vector<float>& x0, std::vector<float>& u, std::vector<float>& v, std::vector<float>& w, float diff, float dt) {

    }

    void add_source(std::vector<float>& x, std::vector<float>& s, float dt) {
    }

    void diffuse(int b, std::vector<float>& x, std::vector<float>& x0, float diff, float dt) {

    }
};

#endif