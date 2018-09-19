#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <cmath>
#include "bvp1d.h"
#include "var_expr.h"

inline int nindex (const Data1D& data, int i, int j, int n)
{
    return i * data.grid.number_of_nodes + data.grid.index(j, n);
}

VarExpr U (const Data1D& data, int i, int j, int n)
{
    VarExpr ve;
    ve.num = 1;
    ve.ind.push_back(nindex(data, i, j, n));
    ve.val.push_back(1);
    ve.rhs = 0;
    return ve;
}

VarExpr DifferenceOperator1D (const Data1D& data, int i, int j, int n)
{
    int M = data.grid.M;
    if (j == 0 && n == 0 || j == M - 1 && n == data.grid.K[j]) {
        // boundary node
        int n1;
        double b0, w0;
        if (n == 0) {
            // left boundary node
            n1 = n + 1;
            b0 = data.b[i][0];
            w0 = data.w[i][0];
        }
        else {  // n == data.grid.K[M - 1]
            // right boundary node
            n1 = n - 1;
            b0 = data.b[i][1];
            w0 = data.w[i][1];
        }
        if (std::isinf(b0)) {
            // Dirichlet BC
            return U(data, i, j, n) - w0;
        }
        else {
            // Neumann or Robin BC
            return data.a[i][j] / data.grid.h[j] *
                    (U(data, i, j, n) - U(data, i, j, n1))
                + b0 * U(data, i, j, n) - w0;
        }
    }
    else if (n == data.grid.K[j]) {
        // rightmost node of the j-th domain
        if (std::isinf(data.G[i][j])) {
            // perfect contact conjugation conditions
            return data.a[i][j] / data.grid.h[j] *
                    (U(data, i, j, n) - U(data, i, j, n - 1))
                + data.a[i][j + 1] / data.grid.h[j + 1] *
                    (U(data, i, j + 1, 0) - U(data, i, j + 1, 1));
        }
        else {
            // imperfect contact conjugation conditions
            return data.a[i][j] / data.grid.h[j] *
                    (U(data, i, j, n) - U(data, i, j, n - 1))
                + data.G[i][j] * (U(data, i, j, n) - U(data, i, j + 1, 0));
        }
    }
    else if (n == 0) {
        // leftmost node of the j-th domain
        if (std::isinf(data.G[i][j - 1])) {
            // perfect contact conjugation conditions
            return U(data, i, j, n) - U(data, i, j - 1, data.grid.K[j - 1]);
        }
        else {
            // imperfect contact conjugation conditions
            return data.a[i][j] / data.grid.h[j]
                    * (U(data, i, j, n) - U(data, i, j, n + 1))
                + data.G[i][j - 1]
                    * (U(data, i, j, n) - U(data, i, j - 1, data.grid.K[j - 1]));
        }
    }
    else {
        // internal node
        return - data.a[i][j] / pow(data.grid.h[j], 2)
            * (U(data, i, j, n - 1) - 2 * U(data, i, j, n)
               + U(data, i, j, n + 1));
    }
}

double DifferenceOperatorValue1D (const Data1D& data,
    const std::vector<GridFunction1D>& u, int i, int j, int n)
{
    int M = data.grid.M;
    if (j == 0 && n == 0 || j == M - 1 && n == data.grid.K[j]) {
        // boundary node
        int n1;
        double b0, w0;
        if (n == 0) {
            // left boundary node
            n1 = n + 1;
            b0 = data.b[i][0];
            w0 = data.w[i][0];
        }
        else {  // n == data.grid.K[M - 1]
            // right boundary node
            n1 = n - 1;
            b0 = data.b[i][1];
            w0 = data.w[i][1];
        }
        if (std::isinf(b0)) {
            // Dirichlet BC
            return 0.0;
        }
        else {
            // Neumann or Robin BC
            return data.a[i][j] / data.grid.h[j] * (u[i](j, n) - u[i](j, n1))
                + b0 * u[i](j, n) - w0;
        }
    }
    else if (n == data.grid.K[j]) {
        // rightmost node of the j-th domain
        // TODO
        throw;
    }
    else if (n == 0) {
        // leftmost node of the j-th domain
        // TODO
        throw;
    }
    else {
        // internal node
        return - data.a[i][j] / pow(data.grid.h[j], 2)
            * (u[i](j, n - 1) - 2 * u[i](j, n) + u[i](j, n + 1));
    }
}

VarExpr get_nonlinear_term (const Data1D& data, int i, int j, int n,
    const mtl::dense_vector<double>& x_old)
{
    VarExpr ve = data.c[i] * U(data, i, j, n) - data.g[i](j, n);
    for (int k = 0; k < data.N; ++k) {
        int ind = nindex(data, k, j, n);
        ve += data.f[i][j][k](x_old[ind]) +
            data.df[i][j][k](x_old[ind]) * (U(data, k, j, n) - x_old[ind]);
    }
    return ve;
}

double get_nonlinear_term_value (const Data1D& data,
    const std::vector<GridFunction1D>& u, int i, int j, int n)
{
    double ans = data.c[i] * u[i](j, n) - data.g[i](j, n);
    for (int k = 0; k < data.N; ++k)
        ans += data.f[i][j][k](u[k](j, n));
    return ans;
}

VarExpr NonlinearOperator1D (const Data1D& data, int i, int j, int n,
    const mtl::dense_vector<double>& x_old)
{
    int M = data.grid.M;
    if (j == 0 && n == 0 || j == M - 1 && n == data.grid.K[j]) {
        // boundary node
        double b0;
        if (n == 0) {
            // left boundary node
            b0 = data.b[i][0];
        }
        else {  // n == data.grid.K[M - 1]
            // right boundary node
            b0 = data.b[i][1];
        }
        if (std::isinf(b0)) {
            // Dirichlet BC
            return 0;
        }
        else {
            // Neumann or Robin BC
            return data.grid.h[j] / 2
                * get_nonlinear_term(data, i, j, n, x_old);
        }
    }
    else if (n == data.grid.K[j]) {
        // rightmost node of the j-th domain
        if (std::isinf(data.G[i][j])) {
            // perfect contact conjugation conditions
            return data.grid.h[j] / 2
                    * get_nonlinear_term(data, i, j, n, x_old)
                + data.grid.h[j + 1] / 2
                    * get_nonlinear_term(data, i, j + 1, 0, x_old);
        }
        else {
            // imperfect contact conjugation conditions
            return data.grid.h[j] / 2
                * get_nonlinear_term(data, i, j, n, x_old);
        }
    }
    else if (n == 0) {
        // leftmost node of the j-th domain
        if (std::isinf(data.G[i][j - 1])) {
            // perfect contact conjugation conditions
            return 0;
        }
        else {
            // imperfect contact conjugation conditions
            return data.grid.h[j] / 2
                * get_nonlinear_term(data, i, j, n, x_old);
        }
    }
    else {
        // internal node
        return get_nonlinear_term(data, i, j, n, x_old);
    }
}

double NonlinearOperatorValue1D (const Data1D& data,
    const std::vector<GridFunction1D>& u, int i, int j, int n)
{
    int M = data.grid.M;
    if (j == 0 && n == 0 || j == M - 1 && n == data.grid.K[j]) {
        // boundary node
        double b0;
        if (n == 0) {
            // left boundary node
            b0 = data.b[i][0];
        }
        else {  // n == data.grid.K[M - 1]
            // right boundary node
            b0 = data.b[i][1];
        }
        if (std::isinf(b0)) {
            // Dirichlet BC
            return 0.0;
        }
        else {
            // Neumann or Robin BC
            return data.grid.h[j] / 2
                * get_nonlinear_term_value(data, u, i, j, n);
        }
    }
    else if (n == data.grid.K[j]) {
        // rightmost node of the j-th domain
        // TODO
        throw;
    }
    else if (n == 0) {
        // leftmost node of the j-th domain
        // TODO
        throw;
    }
    else {
        // internal node
        return get_nonlinear_term_value(data, u, i, j, n);
    }
}

void SolveBVP1D (const Data1D& data, const Parameters1D& param,
    std::vector<GridFunction1D>& sol)
{
    int unknowns = data.N * data.grid.number_of_nodes;
    mtl::compressed2D<double> A(unknowns, unknowns);
    mtl::dense_vector<double> b(unknowns), x(unknowns), x_old(unknowns);
    int elements_per_row = 2 + data.N;  // parameter of the inserter
    // Set the initial guess
    for (int i = 0; i < data.N; ++i) {
        for (int j = 0; j < data.grid.M; ++j) {
            for (int n = 0; n <= data.grid.K[j]; ++n)
                x[nindex(data, i, j, n)] = sol[i](j, n);
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
                for (int j = 0; j < data.grid.M; ++j) {
                    for (int n = 0; n <= data.grid.K[j]; ++n) {
                        VarExpr ve = DifferenceOperator1D(data, i, j, n)
                            + NonlinearOperator1D(data, i, j, n, x_old);
                        int row = nindex(data, i, j, n);
                        b[row] = ve.rhs;
                        for (int s = 0; s < ve.num; ++s)
                            ins(row, ve.ind[s]) << ve.val[s];
                    }
                }
            }
        }
        // Solve the linear system
        if (param.sol_method == SOL_METHOD_MTL) {
            itl::pc::ilu_0<mtl::compressed2D<double> > P(A);
            itl::basic_iteration<double> iter(b, param.max_linear_sys_iterations,
                param.linear_sys_tol);
            // x is equal to the previous guess
            bicgstab(A, x, b, P, iter);
        }
        else {  // sol_method == SOL_METHOD_UMFPACK
            #ifdef MTL_HAS_UMFPACK
                umfpack_solve(A, x, b);
            #else
                std::cerr << "UMFPACK library is needed!" << std::endl;
                throw;
            #endif // MTL_HAS_UMFPACK
        }
        ++num_iterations;
    } while (!(mtl::infinity_norm(x - x_old) < param.Newton_tol ||
        num_iterations >= param.max_Newton_iterations));
    // Set the solution
    for (int i = 0; i < data.N; ++i) {
        for (int j = 0; j < data.grid.M; ++j) {
            for (int n = 0; n <= data.grid.K[j]; ++n)
                sol[i](j, n) = x[nindex(data, i, j, n)];
        }
    }
}

double OperatorValue1D (const Data1D& data,
    const std::vector<GridFunction1D>& u, int i, int j, int n)
{
    return DifferenceOperatorValue1D(data, u, i, j, n)
        + NonlinearOperatorValue1D(data, u, i, j, n);
}
