#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <cmath>
#include "bvp2d.h"
#include "var_expr.h"

inline int nindex (const Data2D& data, int i, int jX, int jY, int nX, int nY)
{
    return i * data.grid.total_number_of_nodes + data.grid.index(jX, jY, nX, nY);
}

VarExpr U(const Data2D& data, int i, int jX, int jY, int nX, int nY)
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
            // left boundary node
            n1[var] = n[var] + 1;
            b0 = data.b[i](var, 0, j[var1], n[var1]);
            w0 = data.w[i](var, 0, j[var1], n[var1]);
        }
        else {  // n[var] == data.grid.K[var][j[var]]
            // right boundary node
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
                + data.G[i][var][jX][jY] * (U(data, i, jX, jY, nX, nY)
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
                + data.G[i][var][j2[VAR_X]][j2[VAR_Y]] * (U(data, i, jX, jY, nX, nY)
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

VarExpr DifferenceOperator2D (const Data2D& data, int i, int jX, int jY, int nX, int nY)
{
    // TODO: Dirichlet BC
    return DifferenceOperator2D_1D(data, VAR_X, i, jX, jY, nX, nY)
        + DifferenceOperator2D_1D(data, VAR_Y, i, jX, jY, nX, nY);
}

VarExpr NonlinearOperator2D (const Data2D& data, int i, int jX, int jY,
    int nX, int nY, const mtl::dense_vector<double>& x_old)
{
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
    mtl::compressed2D<double> A(unknowns, unknowns);
    mtl::dense_vector<double> b(unknowns), x(unknowns), x_old(unknowns);
    int elements_per_row = 4 + data.N;  // parameter of the inserter
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
    do {
        A = 0;
        x_old = x;
        {  // additional block for applying the inserter
            mtl::mat::inserter<mtl::compressed2D<double>,
                mtl::update_plus<double> > ins(A, elements_per_row);
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
                                    ins(row, ve.ind[s]) << ve.val[s];
                            }
                        }
                    }
                }
            }
        }
        // Solve the linear system
        itl::pc::ilu_0<mtl::compressed2D<double> > P(A);
        itl::basic_iteration<double> iter(b, param.max_linear_sys_iterations,
            param.linear_sys_tol);
        // x is equal to the previous guess
        bicgstab(A, x, b, P, iter);
        ++num_iterations;
    } while (!(mtl::infinity_norm(x - x_old) < param.Newton_tol ||
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
