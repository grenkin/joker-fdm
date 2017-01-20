#ifndef BVP2D_H_INCLUDED
#define BVP2D_H_INCLUDED

#include <vector>
#include <array>
#include <cstring>
#include <cstddef>

enum var_t { VAR_X, VAR_Y };

class BoundaryGridFunction2D;

class Grid2D {
    friend BoundaryGridFunction2D;
    std::vector<int> first_index[2];
    std::vector<double> start_coord[2];
    int number_of_nodes[2];
public:
    int M[2];
    std::array<std::vector<double>, 2> L, h;
    std::array<std::vector<int>, 2> K;
    int total_number_of_nodes;
    Grid2D (const std::array<std::vector<double>, 2>& _L,
        const std::array<std::vector<int>, 2>& _K)
        : L(_L), K(_K)
    {
        for (int d = 0; d < 2; ++d) {
            M[d] = L[d].size();
            h[d].resize(M[d]);
            for (int j = 0; j < M[d]; ++j)
                h[d][j] = L[d][j] / K[d][j];
            start_coord[d].resize(M[d]);
            start_coord[d][0] = 0;
            for (int j = 1; j < M[d]; ++j)
                start_coord[d][j] = start_coord[d][j - 1] + L[d][j - 1];
            first_index[d].resize(M[d]);
            first_index[d][0] = 0;
            for (int j = 1; j < M[d]; ++j)
                first_index[d][j] = first_index[d][j - 1] + K[d][j - 1] + 1;
            number_of_nodes[d] = first_index[d][M[d] - 1] + K[d][M[d] - 1] + 1;
        }
        total_number_of_nodes = number_of_nodes[VAR_X] * number_of_nodes[VAR_Y];
    }
    int index (int jX, int jY, int nX, int nY) const
    {
        return (first_index[VAR_X][jX] + nX) * number_of_nodes[VAR_Y]
            + first_index[VAR_Y][jY] + nY;
    }
    double coordX (int jX, int nX) const
    {
        return start_coord[VAR_X][jX] + h[VAR_X][jX] * nX;
    }
    double coordY (int jY, int nY) const
    {
        return start_coord[VAR_Y][jY] + h[VAR_Y][jY] * nY;
    }
};

class GridFunction2D {
    Grid2D grid;
    double* v;
public:
    GridFunction2D ()
        : grid({std::vector<double>(1, 0.0), std::vector<double>(1, 0.0)},
               {std::vector<int>(1, 0), std::vector<int>(1, 0)})
    {
        v = NULL;
    }
    GridFunction2D (const Grid2D& _grid)
        : grid(_grid)
    {
		v = new double[grid.total_number_of_nodes];
    }
    GridFunction2D (const GridFunction2D& gf)
        : grid(gf.grid)
    {
        v = new double[grid.total_number_of_nodes];
        memcpy(v, gf.v, grid.total_number_of_nodes * sizeof(double));
    }
    ~GridFunction2D ()
    {
        delete[] v;
    }
    void set_grid (const Grid2D& _grid)
    {
        grid = _grid;
        v = new double[grid.total_number_of_nodes];
    }
    double& operator() (int jX, int jY, int nX, int nY)
    {
        return v[grid.index(jX, jY, nX, nY)];
    }
    const double& operator() (int jX, int jY, int nX, int nY) const
    {
        return v[grid.index(jX, jY, nX, nY)];
    }
};

class BoundaryGridFunction2D {
    Grid2D grid;
    double* v;
    int number_of_elements()
    {
        return 2 * (grid.number_of_nodes[VAR_X] + grid.number_of_nodes[VAR_Y]);
    }
    int index(var_t var, int s, int j, int n) const
    {
        if (var == VAR_X) {
            return s * grid.number_of_nodes[VAR_Y]
                + grid.first_index[VAR_Y][j] + n;
        }
        else {  // var == VAR_Y
            return 2 * grid.number_of_nodes[VAR_Y]
                + s * grid.number_of_nodes[VAR_X] + grid.first_index[VAR_X][j] + n;
        }
    }
public:
    BoundaryGridFunction2D ()
        : grid({std::vector<double>(1, 0.0), std::vector<double>(1, 0.0)},
               {std::vector<int>(1, 0), std::vector<int>(1, 0)})
    {
        v = NULL;
    }
    BoundaryGridFunction2D (const Grid2D& _grid)
        : grid(_grid)
    {
		v = new double[number_of_elements()];
    }
    BoundaryGridFunction2D (const BoundaryGridFunction2D& gf)
        : grid(gf.grid)
    {
        v = new double[number_of_elements()];
        memcpy(v, gf.v, number_of_elements() * sizeof(double));
    }
    ~BoundaryGridFunction2D ()
    {
        delete[] v;
    }
    void set_grid (const Grid2D& _grid)
    {
        grid = _grid;
        v = new double[number_of_elements()];
    }
    double& operator() (var_t var, int s, int j, int n)
    {
        return v[index(var, s, j, n)];
    }
    const double& operator() (var_t var, int s, int j, int n) const
    {
        return v[index(var, s, j, n)];
    }
};

typedef double (*NonlinearFunction) (double);

struct Data2D {
    Grid2D grid;
    int N;  // number of equations
    std::vector<std::vector<std::vector<double> > > a;  // diffusion coefficients
    std::vector<std::array<std::vector<std::vector<double> >, 2> > G;  // coefficients in conjugation conditions
    std::vector<BoundaryGridFunction2D> b, w;  // coefficients in boundary conditions
    std::vector<std::vector<std::vector<std::vector<NonlinearFunction> > > > f, df;  // nonlinear terms functions
    std::vector<GridFunction2D> g;  // right-hand sides
    Data2D (int _N, const Grid2D& _grid)
        : grid(_grid), N(_N)
    {
        a.resize(N);
        G.resize(N);
        b.resize(N);
        w.resize(N);
        for (int i = 0; i < N; ++i) {
            a[i].resize(grid.M[VAR_X]);
            for (int jX = 0; jX < grid.M[VAR_X]; ++jX)
                a[i][jX].resize(grid.M[VAR_Y]);
            G[i][VAR_X].resize(grid.M[VAR_X] - 1);
            for (int jX = 0; jX < grid.M[VAR_X] - 1; ++jX)
                G[i][VAR_X][jX].resize(grid.M[VAR_Y]);
            G[i][VAR_Y].resize(grid.M[VAR_X]);
            for (int jX = 0; jX < grid.M[VAR_X]; ++jX)
                G[i][VAR_Y][jX].resize(grid.M[VAR_Y] - 1);
            b[i].set_grid(grid);
            w[i].set_grid(grid);
        }
        f.resize(N);
        for (int i = 0; i < N; ++i) {
            f[i].resize(grid.M[VAR_X]);
            for (int jX = 0; jX < grid.M[VAR_X]; ++jX) {
                f[i][jX].resize(grid.M[VAR_Y]);
                for (int jY = 0; jY < grid.M[VAR_Y]; ++jY)
                    f[i][jX][jY] = std::vector<NonlinearFunction>(N, NULL);
            }
        }
        df = f;
        g.resize(N);
        for (int i = 0; i < N; ++i)
            g[i].set_grid(grid);
    }
};

struct Parameters2D {
    double Newton_tol, linear_sys_tol;
    int max_Newton_iterations, max_linear_sys_iterations;
    Parameters2D ()
    {
        Newton_tol = 1e-8;  // absolute error in Newton's method
        max_Newton_iterations = 100;  // maximal number of iterations of Newton's method
        linear_sys_tol = 1e-7;  // relative error in the linear system solver
        max_linear_sys_iterations = 10000;  // maximal number of iterations in the linear system solver
    }
};

void SolveBVP2D (const Data2D&, const Parameters2D&,
    std::vector<GridFunction2D>&);

#endif // BVP2D_H_INCLUDED
