#include <mtl/matrix.h>
#include <itl/preconditioner/ilu.h>
#include <itl/interface/mtl.h>
#include <itl/krylov/bicgstab.h>
#include <cmath>
#include "bvp2d.h"
#include "var_expr.h"

inline int nindex (const Data2D& data, int i, int jX, int jY, int nX, int nY)
{
    return i * data.grid.total_number_of_nodes + data.grid.index(jX, jY, nX, nY);
}

VarExpr U (const Data2D& data, int i, int jX, int jY, int nX, int nY)
{
    VarExpr ve;
    ve.num = 1;
    ve.ind.push_back(nindex(data, i, jX, jY, nX, nY));
    ve.val.push_back(1);
    ve.rhs = 0;
    return ve;
}

VarExpr DifferenceOperator2D_1D (const Data2D& data, var_t var,
    int i, int jX, int jY, int nX, int nY)
    // except the case of Dirichlet BC
{
    std::array<int, 2> j = {jX, jY}, n = {nX, nY};
    if (j[var] == 0 && n[var] == 0 || j[var] == data.grid.M[var] - 1
        && n[var] == data.grid.K[var][j[var]])
    {
        // boundary node
        var_t var1 = (var_t)(1 - var);
        std::array<int, 2> n1 = n;
        double b0, w0;
        if (n[var] == 0) {
            // left or bottom boundary node
            n1[var] = n[var] + 1;
            b0 = data.b[i](var, 0, j[var1], n[var1]);
            w0 = data.w[i](var, 0, j[var1], n[var1]);
        }
        else {  // n[var] == data.grid.K[var][j[var]]
            // right or top boundary node
            n1[var] = n[var] - 1;
            b0 = data.b[i](var, 1, j[var1], n[var1]);
            w0 = data.w[i](var, 1, j[var1], n[var1]);
        }
        if (std::isinf(b0)) {
            // Dirichlet BC
            throw;
        }
        else {
            return 2 * data.a[i][jX][jY] / pow(data.grid.h[var][j[var]], 2)
                * (U(data, i, jX, jY, nX, nY) - U(data, i, jX, jY, n1[VAR_X], n1[VAR_Y]))
                + 2 / data.grid.h[var][j[var]] * (b0 * U(data, i, jX, jY, nX, nY) - w0);
        }
    }
    else if (n[var] == data.grid.K[var][j[var]]) {
        // rightmost or topmost node of the j-th domain
        if (std::isinf(data.G[i][var][jX][jY])) {
            // perfect contact conjugation conditions
            // TODO; choose a large G instead of infinity
            throw;
        }
        else {
            // imperfect contact conjugation conditions
            std::array<int, 2> n1 = n, j2 = j, n2 = n;
            n1[var] = n[var] - 1;
            j2[var] = j[var] + 1;
            n2[var] = 0;
            return 2 * data.a[i][jX][jY] / pow(data.grid.h[var][j[var]], 2)
                * (U(data, i, jX, jY, nX, nY) - U(data, i, jX, jY, n1[VAR_X], n1[VAR_Y]))
                + 2 * data.G[i][var][jX][jY] / data.grid.h[var][j[var]]
                    * (U(data, i, jX, jY, nX, nY)
                        - U(data, i, j2[VAR_X], j2[VAR_Y], n2[VAR_X], n2[VAR_Y]));
        }
    }
    else if (n[var] == 0) {
        // leftmost or bottommost node of the j-th domain
        std::array<int, 2> j2 = j;
        j2[var] = j[var] - 1;
        if (std::isinf(data.G[i][var][j2[VAR_X]][j2[VAR_Y]])) {
            // perfect contact conjugation conditions
            // TODO; choose a large G instead of infinity
            throw;
        }
        else {
            // imperfect contact conjugation conditions
            std::array<int, 2> n1 = n, n2 = n;
            n1[var] = n[var] + 1;
            n2[var] = data.grid.K[var][j2[var]];
            return 2 * data.a[i][jX][jY] / pow(data.grid.h[var][j[var]], 2)
                * (U(data, i, jX, jY, nX, nY) - U(data, i, jX, jY, n1[VAR_X], n1[VAR_Y]))
                + 2 * data.G[i][var][j2[VAR_X]][j2[VAR_Y]] / data.grid.h[var][j[var]]
                    * (U(data, i, jX, jY, nX, nY)
                        - U(data, i, j2[VAR_X], j2[VAR_Y], n2[VAR_X], n2[VAR_Y]));
        }
    }
    else {
        // internal node
        std::array<int, 2> n1 = n, n2 = n;
        n1[var] = n[var] - 1;
        n2[var] = n[var] + 1;
        return - data.a[i][jX][jY] / pow(data.grid.h[var][j[var]], 2)
            * (U(data, i, jX, jY, n1[VAR_X], n1[VAR_Y])
               - 2 * U(data, i, jX, jY, nX, nY)
               + U(data, i, jX, jY, n2[VAR_X], n2[VAR_Y]));
    }
}

bool equal(double a, double b)
{
    return abs(a - b) < 1e-12;
}

bool is_Dirichlet (const Data2D& data, int i, int jX, int jY, int nX, int nY, double& w0)
{
    // Note: equal nodes in conjugation conditions are not checked
    std::array<int, 2> j = {jX, jY}, n = {nX, nY};
    bool w0_set = false;
    for (int var_index = 0; var_index < 2; ++var_index) {
        var_t var = (var_t)var_index;
        if (j[var] == 0 && n[var] == 0 || j[var] == data.grid.M[var] - 1
            && n[var] == data.grid.K[var][j[var]])
        {
            var_t var1 = (var_t)(1 - var_index);
            double b1, w1;
            if (n[var] == 0) {
                // left or bottom boundary node
                b1 = data.b[i](var, 0, j[var1], n[var1]);
                w1 = data.w[i](var, 0, j[var1], n[var1]);
            }
            else {  // n[var] == data.grid.K[var][j[var]]
                // right or top boundary node
                b1 = data.b[i](var, 1, j[var1], n[var1]);
                w1 = data.w[i](var, 1, j[var1], n[var1]);
            }
            if (std::isinf(b1)) {
                if (w0_set) {
                    if (!equal(w0, w1))
                        throw ENotConsistentDirichletBC(i, jX, jY, nX, nY);
                }
                else {
                    w0 = w1;
                    w0_set = true;
                }
            }
        }
    }
    return w0_set;
}

VarExpr DifferenceOperator2D (const Data2D& data, int i, int jX, int jY, int nX, int nY)
{
    double w0;
    if (is_Dirichlet(data, i, jX, jY, nX, nY, w0))
        return U(data, i, jX, jY, nX, nY) - w0;
    else {
        return DifferenceOperator2D_1D(data, VAR_X, i, jX, jY, nX, nY)
            + DifferenceOperator2D_1D(data, VAR_Y, i, jX, jY, nX, nY);
    }
}

VarExpr NonlinearOperator2D (const Data2D& data, int i, int jX, int jY,
    int nX, int nY, const mtl::dense1D<double>& x_old)
{
    double w0;
    if (is_Dirichlet(data, i, jX, jY, nX, nY, w0))
        return 0;
    VarExpr ve = - data.g[i](jX, jY, nX, nY);
    for (int k = 0; k < data.N; ++k) {
        int ind = nindex(data, k, jX, jY, nX, nY);
        ve += data.f[i][jX][jY][k](x_old[ind]) +
            data.df[i][jX][jY][k](x_old[ind]) * (U(data, k, jX, jY, nX, nY) - x_old[ind]);
    }
    return ve;
}

void SolveBVP2D (const Data2D& data, const Parameters2D& param,
    std::vector<GridFunction2D>& sol)
{
    int unknowns = data.N * data.grid.total_number_of_nodes;
    typedef mtl::matrix<double, mtl::rectangle<>,
        mtl::array<mtl::compressed<> >, mtl::row_major>::type Matrix;
    Matrix A(unknowns, unknowns);
    mtl::dense1D<double> b(unknowns), x(unknowns), x_old(unknowns);
    // Set the initial guess
    for (int i = 0; i < data.N; ++i) {
        for (int jX = 0; jX < data.grid.M[VAR_X]; ++jX) {
            for (int nX = 0; nX <= data.grid.K[VAR_X][jX]; ++nX) {
                for (int jY = 0; jY < data.grid.M[VAR_Y]; ++jY) {
                    for (int nY = 0; nY <= data.grid.K[VAR_Y][jY]; ++nY)
                        x[nindex(data, i, jX, jY, nX, nY)] = sol[i](jX, jY, nX, nY);
                }
            }
        }
    }
    // Apply Newton's method
    int num_iterations = 0;
    double max_diff;
    do {
        for (int i = 0; i < unknowns; ++i)
            x_old[i] = x[i];
        for (int i = 0; i < data.N; ++i) {
            for (int jX = 0; jX < data.grid.M[VAR_X]; ++jX) {
                for (int nX = 0; nX <= data.grid.K[VAR_X][jX]; ++nX) {
                    for (int jY = 0; jY < data.grid.M[VAR_Y]; ++jY) {
                        for (int nY = 0; nY <= data.grid.K[VAR_Y][jY]; ++nY) {
                            VarExpr ve = DifferenceOperator2D(data, i, jX, jY, nX, nY)
                                + NonlinearOperator2D(data, i, jX, jY, nX, nY, x_old);
                            int row = nindex(data, i, jX, jY, nX, nY);
                            b[row] = ve.rhs;
                            for (int s = 0; s < ve.num; ++s)
                                A(row, ve.ind[s]) = 0.0;
                            for (int s = 0; s < ve.num; ++s)
                                A(row, ve.ind[s]) += ve.val[s];
                        }
                    }
                }
            }
        }
        // Solve the linear system
        itl::ILU<Matrix> precond(A);
        itl::basic_iteration<double> iter(b, param.max_linear_sys_iterations,
            param.linear_sys_tol);
        // x is equal to the previous guess
        bicgstab(A, x, b, precond(), iter);
        if (iter.error_code())
            throw;
        ++num_iterations;
        max_diff = 0.0;
        for (int i = 0; i < unknowns; ++i)
            max_diff = fmax(max_diff, fabs(x[i] - x_old[i]));
    } while (!(max_diff < param.Newton_tol ||
        num_iterations >= param.max_Newton_iterations));
    // Set the solution
    for (int i = 0; i < data.N; ++i) {
        for (int jX = 0; jX < data.grid.M[VAR_X]; ++jX) {
            for (int nX = 0; nX <= data.grid.K[VAR_X][jX]; ++nX) {
                for (int jY = 0; jY < data.grid.M[VAR_Y]; ++jY) {
                    for (int nY = 0; nY <= data.grid.K[VAR_Y][jY]; ++nY)
                        sol[i](jX, jY, nX, nY) = x[nindex(data, i, jX, jY, nX, nY)];
                }
            }
        }
    }
}
