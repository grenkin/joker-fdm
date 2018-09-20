/*
An initial-boundary value problem in 1D
3u_t - u_xx + 2v = 2,  2v_t - 4v_xx + 10v - 12u = 10,  x \in (0, pi)
-u_x(0, t) + 2u(0, t) = -exp(-t)
u(0, t) = 0,  u(pi, t) = 0
v(0, t) = 1,  v(pi, t) = 1
u(x, 0) = sin(x),  v(x, 0) = sin(x) + 1
Exact solution: u(x, t) = exp(-t)*sin(x), v(x, t) = exp(-t)*sin(x) + 1

Crank-Nicolson scheme is used:
2 u^m / tau + Au^m = g^{m-1} + g^m + 2 u^{m-1} / tau - Au^{m-1}
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp1d.h>

using namespace std;

const double pi = 3.14159265358979323846;

double zero (double x)
{
    return 0;
}

double two_ident (double x)
{
    return 2 * x;
}

double two (double x)
{
    return 2;
}

double ten_ident (double x)
{
    return 10 * x;
}

double ten (double x)
{
    return 10;
}

double minus_twelve_ident (double x)
{
    return -12 * x;
}

double minus_twelve (double x)
{
    return -12;
}

double gfun1 (double x, double t)
{
    return 2;
}

double gfun2 (double x, double t)
{
    return 10;
}

double initfun1 (double x)
{
    return sin(x);
}

double initfun2 (double x)
{
    return sin(x) + 1;
}

double exact1 (double x, double t)
{
    return exp(-t) * sin(x);
}

double exact2 (double x, double t)
{
    return exp(-t) * sin(x) + 1;
}

int main ()
{
    double L0 = pi;
    double T = 10;
    int K0 = 10;
    int tnum = 10;
    double rms_old;
    for (int test = 0; test <= 5; ++test) {
        vector<double> L(1);  L[0] = L0;
        vector<int> K(1);  K[0] = K0;
        Grid1D grid(L, K);
        Data1D data(2, grid);
        data.a[0][0] = 1;
        data.a[1][0] = 4;
        data.b[0][0] = INFINITY;  data.w[0][0] = 0;
        data.b[0][1] = INFINITY;  data.w[0][1] = 0;
        data.b[1][0] = INFINITY;  data.w[1][0] = 1;
        data.b[1][1] = INFINITY;  data.w[1][1] = 1;
        data.f[0][0][0] = data.df[0][0][0] = zero;
        data.f[0][0][1] = two_ident;  data.df[0][0][1] = two;
        data.f[1][0][0] = minus_twelve_ident;  data.df[1][0][0] = minus_twelve;
        data.f[1][0][1] = ten_ident;  data.df[1][0][1] = ten;
        vector<double> c(2);
        c[0] = 3;  c[1] = 2;
        double tau = T / tnum;
        for (int i = 0; i < 2; ++i)
            data.c[i] = 2 * c[i] / tau;
        vector<GridFunction1D> sol(2);
        vector<TimeGridFunction1D> res(2);
        for (int i = 0; i < 2; ++i) {
            sol[i].set_grid(grid);
            res[i].set_grid(grid, tnum);
        }

        // initial conditions
        for (int n = 0; n <= K0; ++n) {
            sol[0](0, n) = initfun1(grid.coord(0, n));
            sol[1](0, n) = initfun2(grid.coord(0, n));
            for (int i = 0; i < 2; ++i)
                res[i](0, 0, n) = sol[i](0, n);
        }

        for (int m = 1; m <= tnum; ++m) {
            double t = m * tau, t_prev = (m - 1) * tau;
            for (int n = 0; n <= K0; ++n) {
                double x = grid.coord(0, n);
                data.g[0](0, n) = gfun1(x, t_prev);
                data.g[1](0, n) = gfun2(x, t_prev);
            }
            for (int i = 0; i < 2; ++i)
                data.c[i] = - data.c[i];
            for (int n = 0; n <= K0; ++n) {
                double x = grid.coord(0, n);
                data.g[0](0, n) = gfun1(x, t)
                    - OperatorValue1D(data, sol, 0, 0, n);
                data.g[1](0, n) = gfun2(x, t)
                    - OperatorValue1D(data, sol, 1, 0, n);
            }
            for (int i = 0; i < 2; ++i)
                data.c[i] = - data.c[i];
            Parameters1D param;
            param.max_Newton_iterations = 1;
            SolveBVP1D(data, param, sol);
            for (int i = 0; i < 2; ++i) {
                for (int n = 0; n <= K0; ++n)
                    res[i](m, 0, n) = sol[i](0, n);
            }
        }

        double rms = 0;
        for (int m = 0; m <= tnum; ++m) {
            for (int n = 0; n <= K0; ++n)
                rms += pow(res[0](m, 0, n) - exact1(grid.coord(0, n), m * tau), 2)
                    + pow(res[1](m, 0, n) - exact2(grid.coord(0, n), m * tau), 2);
        }
        rms = sqrt(rms / 2 / (K0 + 1) / (tnum + 1));

        cout << "h = " << grid.h[0] << endl;
        cout << "tau = " << tau << endl;
        cout << "rms = " << rms << endl;
        if (test > 0)
            cout << "rms_old / rms = " << rms_old / rms << "\n\n";
        rms_old = rms;
        K0 *= 2;
        tnum *= 2;
    }

    return 0;
}
