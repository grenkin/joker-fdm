#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>
#include <cmath>
#include "bvp1d.h"

inline int nindex (const Data1D& data, int i, int j, int n)
{
    return i * data.grid.number_of_nodes + data.grid.index(j, n);
}

void DifferenceOperator1D (const Data1D& data, int i, int j, int n,
    int& vnum, int vind[], double vval[], double& rhs)
{
    int M = data.grid.M;
    if (j == 0 && n == 0 || j == M - 1 && n == data.grid.K[j]) {
        // boundary node
        int n1;
        double b0, v0;
        if (n == 0) {
            // left boundary node
            n1 = n + 1;
            b0 = data.b[i][0];
            v0 = data.v[i][0];
        }
        else {  // n == data.grid.K[M - 1]
            // right boundary node
            n1 = n - 1;
            b0 = data.b[i][1];
            v0 = data.v[i][1];
        }
        if (std::isinf(b0)) {
            // Dirichlet BC
            vnum = 1;
            vind[0] = nindex(data, i, j, n);
            vval[0] = 1;
            rhs = v0;
        }
        else {
            // Neumann or Robin BC
            vnum = 2;
            double c = data.a[i][j] / data.grid.h[j];
            vind[0] = nindex(data, i, j, n);
            vval[0] = c + b0;
            vind[1] = nindex(data, i, j, n1);
            vval[1] = c * -1;
            rhs = b0 * v0;
        }
    }
    else if (n == data.grid.K[j]) {
        // rightmost node of the j-th domain
        if (std::isinf(data.G[i][j])) {
            // perfect contact conjugation conditions
            vnum = 4;
            double c1 = data.a[i][j] / data.grid.h[j];
            double c2 = data.a[i][j + 1] / data.grid.h[j + 1];
            vind[0] = nindex(data, i, j, n);
            vval[0] = c1;
            vind[1] = nindex(data, i, j, n - 1);
            vval[1] = c1 * -1;
            vind[2] = nindex(data, i, j + 1, 0);
            vval[2] = c2;
            vind[3] = nindex(data, i, j + 1, 1);
            vval[3] = c2 * -1;
            rhs = 0;
        }
        else {
            // imperfect contact conjugation conditions
            vnum = 3;
            double c = data.a[i][j] / data.grid.h[j];
            vind[0] = nindex(data, i, j, n);
            vval[0] = c + data.G[i][j];
            vind[1] = nindex(data, i, j, n - 1);
            vval[1] = c * -1;
            vind[2] = nindex(data, i, j + 1, 0);
            vval[2] = data.G[i][j] * -1;
            rhs = 0;
        }
    }
    else if (n == 0) {
        // leftmost node of the j-th domain
        if (std::isinf(data.G[i][j - 1])) {
            // perfect contact conjugation conditions
            vnum = 2;
            vind[0] = nindex(data, i, j, n);
            vval[0] = 1;
            vind[1] = nindex(data, i, j - 1, data.grid.K[j - 1]);
            vval[1] = -1;
            rhs = 0;
        }
        else {
            // imperfect contact conjugation conditions
            vnum = 3;
            double c = data.a[i][j] / data.grid.h[j];
            vind[0] = nindex(data, i, j, n);
            vval[0] = c + data.G[i][j - 1];
            vind[1] = nindex(data, i, j, n + 1);
            vval[1] = c * -1;
            vind[2] = nindex(data, i, j - 1, data.grid.K[j - 1]);
            vval[2] = data.G[i][j - 1] * -1;
            rhs = 0;
        }
    }
    else {
        // internal node
        vnum = 3;
        double c = - data.a[i][j] / pow(data.grid.h[j], 2);
        vind[0] = nindex(data, i, j, n - 1);
        vval[0] = c;
        vind[1] = nindex(data, i, j, n);
        vval[1] = c * -2;
        vind[2] = nindex(data, i, j, n + 1);
        vval[2] = c;
        rhs = 0;
    }
}

void get_nonlinear_term (const Data1D& data, int i, int j, int n,
    const mtl::dense_vector<double>& x_old,
    int vind[], double vval[], double& rhs)
{
    rhs = data.g[i](j, n);
    for (int k = 0; k < data.N; ++k) {
        vind[k] = nindex(data, k, j, n);
        vval[k] = data.df[i][j][k](x_old[vind[k]]);
        rhs -= data.f[i][j][k](x_old[vind[k]])
            - data.df[i][j][k](x_old[vind[k]]) * x_old[vind[k]];
    }
}

void NonlinearOperator1D (const Data1D& data, int i, int j, int n,
    const mtl::dense_vector<double>& x_old,
    int& vnum, int vind[], double vval[], double& rhs)
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
            vnum = 0;
            rhs = 0;
        }
        else {
            // Neumann or Robin BC
            vnum = data.N;
            get_nonlinear_term(data, i, j, n, x_old, vind, vval, rhs);
            for (int k = 0; k < data.N; ++k)
                vval[k] *= data.grid.h[j] / 2;
            rhs *= data.grid.h[j] / 2;
        }
    }
    else if (n == data.grid.K[j]) {
        // rightmost node of the j-th domain
        if (std::isinf(data.G[i][j])) {
            // perfect contact conjugation conditions
            vnum = data.N;
            double *vval1 = new double[data.N], *vval2 = new double[data.N];
            double rhs1, rhs2;
            get_nonlinear_term(data, i, j, n, x_old, vind, vval1, rhs1);
            get_nonlinear_term(data, i, j + 1, 0, x_old, vind, vval2, rhs2);
            for (int k = 0; k < data.N; ++k) {
                vval[k] = vval1[k] * data.grid.h[j] / 2
                    + vval2[k] * data.grid.h[j + 1] / 2;
            }
            rhs = rhs1 * data.grid.h[j] / 2 + rhs2 * data.grid.h[j + 1] / 2;
            delete[] vval1;
            delete[] vval2;
        }
        else {
            // imperfect contact conjugation conditions
            vnum = data.N;
            get_nonlinear_term(data, i, j, n, x_old, vind, vval, rhs);
            for (int k = 0; k < data.N; ++k)
                vval[k] *= data.grid.h[j] / 2;
            rhs *= data.grid.h[j] / 2;
        }
    }
    else if (n == 0) {
        // leftmost node of the j-th domain
        if (std::isinf(data.G[i][j - 1])) {
            // perfect contact conjugation conditions
            vnum = 0;
            rhs = 0;
        }
        else {
            // imperfect contact conjugation conditions
            vnum = data.N;
            get_nonlinear_term(data, i, j, n, x_old, vind, vval, rhs);
            for (int k = 0; k < data.N; ++k)
                vval[k] *= data.grid.h[j] / 2;
            rhs *= data.grid.h[j] / 2;
        }
    }
    else {
        // internal node
        vnum = data.N;
        get_nonlinear_term(data, i, j, n, x_old, vind, vval, rhs);
    }
}

void SolveBVP1D (const Data1D& data, const Parameters& param,
    std::vector<GridFunction1D>& sol)
{
    int unknowns = data.N * data.grid.number_of_nodes;
    mtl::compressed2D<double> A(unknowns, unknowns);
    mtl::dense_vector<double> b(unknowns), x(unknowns), x_old(unknowns);
    int elements_per_row = 3 + data.N;  // parameter of the inserter
    int max_elements_per_row = 4 + data.N;
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
            int vnum;
            int vind[max_elements_per_row];
            double vval[max_elements_per_row];
            double rhs;
            for (int i = 0; i < data.N; ++i) {
                for (int j = 0; j < data.grid.M; ++j) {
                    for (int n = 0; n <= data.grid.K[j]; ++n) {
                        int row = nindex(data, i, j, n);
                        b[row] = 0;
                        DifferenceOperator1D(data, i, j, n,
                            vnum, vind, vval, rhs);
                        for (int s = 0; s < vnum; ++s)
                            ins(row, vind[s]) << vval[s];
                        b[row] += rhs;
                        NonlinearOperator1D(data, i, j, n, x_old,
                            vnum, vind, vval, rhs);
                        for (int s = 0; s < vnum; ++s)
                            ins(row, vind[s]) << vval[s];
                        b[row] += rhs;
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
        for (int j = 0; j < data.grid.M; ++j) {
            for (int n = 0; n <= data.grid.K[j]; ++n)
                sol[i](j, n) = x[nindex(data, i, j, n)];
        }
    }
}
