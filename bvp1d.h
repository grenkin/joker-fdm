#ifndef BVP_1D_H_INCLUDED
#define BVP_1D_H_INCLUDED

#include <vector>
#include <array>
#include <cstring>
#include <cstddef>

class Grid1D {
    std::vector<int> first_index;  // index of the first node in a subdomain
    std::vector<double> left_coord;  // the leftmost coordinate of a subdomain
public:
    int M;  // number of subdomains
    std::vector<double> L;  // lengths of subdomains
    std::vector<int> K;  // numbers of grid intervals in each subdomain
    std::vector<double> h;  // grid steps
    int number_of_nodes;
    Grid1D (const std::vector<double>& _L, const std::vector<int>& _K)
        : M(_L.size()), L(_L), K(_K)
    {
        h.resize(M);
        for (int j = 0; j < M; ++j)
            h[j] = L[j] / K[j];
        first_index.resize(M);
        first_index[0] = 0;
        for (int j = 1; j < M; ++j)
            first_index[j] = first_index[j - 1] + K[j - 1] + 1;
        left_coord.resize(M);
        left_coord[0] = 0;
        for (int j = 1; j < M; ++j)
            left_coord[j] = left_coord[j - 1] + L[j - 1];
        number_of_nodes = first_index[M - 1] + K[M - 1] + 1;
    }
    int index (int j, int n) const
    {
        return first_index[j] + n;
    }
    double coord (int j, int n) const
    {
        return left_coord[j] + h[j] * n;
    }
};

class GridFunction1D {
    Grid1D grid;
    double *v;
public:
    GridFunction1D ()
        : grid(std::vector<double>(1, 0.0), std::vector<int>(1, 0))
    {
        v = NULL;
    }
    GridFunction1D (const GridFunction1D& gf)
        : grid(gf.grid)
    {
        v = new double[grid.number_of_nodes];
        memcpy(v, gf.v, grid.number_of_nodes * sizeof(double));
    }
    ~GridFunction1D ()
    {
        delete[] v;
    }
    void set_grid (const Grid1D& _grid)
    {
        grid = _grid;
        v = new double[grid.number_of_nodes];
    }
    double& operator() (int j, int n)
    {
        return v[grid.index(j, n)];
    }
    const double& operator() (int j, int n) const
    {
        return v[grid.index(j, n)];
    }
};

typedef double (*NonlinearFunction) (double);

struct Data1D {
    Grid1D grid;
    int N;  // number of equations
    std::vector<std::vector<double> > a;  // diffusion coefficients
    std::vector<std::vector<double> > G;  // coefficients in conjugation conditions
    std::vector<std::array<double, 2> > b, v;  // coefficients in boundary conditions
    std::vector<std::vector<std::vector<NonlinearFunction> > > f, df;  // nonlinear terms functions
    std::vector<GridFunction1D> g;  // right-hand sides
    Data1D (int _N, const Grid1D& _grid)
        : grid(_grid), N(_N)
    {
        int M = grid.M;
        a.resize(N);
        G.resize(N);
        for (int i = 0; i < N; ++i) {
            a[i].resize(M);
            G[i].resize(M - 1);
            b.resize(N);
            v.resize(N);
        }
        f.resize(N);
        for (int i = 0; i < N; ++i) {
            f[i].resize(M);
            for (int j = 0; j < M; ++j)
                f[i][j] = std::vector<NonlinearFunction>(N, NULL);
        }
        df = f;
        g.resize(N);
        for (int i = 0; i < N; ++i)
            g[i].set_grid(grid);
    }
};

struct Parameters {
    double Newton_tol, linear_sys_tol;
    int max_Newton_iterations, max_linear_sys_iterations;
    Parameters ()
    {
        Newton_tol = 1e-8;  // absolute error in Newton's method
        max_Newton_iterations = 100;  // maximal number of iterations of Newton's method
        linear_sys_tol = 1e-7;  // relative error in the linear system solver
        max_linear_sys_iterations = 10000;  // maximal number of iterations in the linear system solver
    }
};

void SolveBVP1D (const Data1D&, const Parameters&,
    std::vector<GridFunction1D>&);

#endif // BVP_1D_H_INCLUDED
