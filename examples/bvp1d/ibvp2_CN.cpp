/*
An initial-boundary value problem in 1D
u_t - u_xx = 2,  -v_xx - u = 0,  x \in (0, 1)
-u_x(0, t) + u(0, t) = 2t
u_x(1, t) + u(1, t) = 2t
-v_x(0, t) + v(0, t) = 0
v_x(1, t) + v(1, t) = 0
u(x, 0) = 0,  v(x, 0) = 0
Exact solution: u(x, t) = 2t, v(x, t) = t*(1 + x - x^2)

Crank-Nicolson scheme is used.
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp1d.h>

using namespace std;

double zero (double x)
{
    return 0;
}

double minus_ident (double x)
{
    return - x;
}

double minus_one (double x)
{
    return -1;
}

double exact1 (double x, double t)
{
    return 2 * t;
}

double exact2 (double x, double t)
{
    return t * (1 + x - x * x);
}

double wfun (double t)
{
    return 2 * t;
}

double gfun0 (double x, double t)
{
    return 2.0;
}

double gfun1 (double x, double t)
{
    return 0.0;
}

double (*gfun[2])(double, double) = {gfun0, gfun1};

int main() {
    double L0 = 1;
    double T = 10;
    int K0 = 10;
    int tnum = 10;
    double rms_old;
    for (int test = 0; test <= 5; ++test) {
        vector<double> L(1);  L[0] = L0;
        vector<int> K(1);  K[0] = K0;
        Grid1D grid(L, K);
        Data1D data(2, grid);
        data.a[0][0] = data.a[1][0] = 1;
        data.b[0][0] = 1;  // data.w[0][0] will be defined later
        data.b[0][1] = 1;  // data.w[0][1] will be defined later
        data.b[1][0] = data.b[1][1] = 1;
        data.w[1][0] = data.w[1][1] = 0;
        data.f[0][0][0] = data.df[0][0][0] = zero;
        data.f[0][0][1] = data.df[0][0][1] = zero;
        data.f[1][0][0] = minus_ident;
        data.df[1][0][0] = minus_one;
        data.f[1][0][1] = data.df[1][0][1] = zero;
        vector<double> c(2);
        c[0] = 1; c[1] = 0;
        double tau = T / tnum;
        for (int i = 0; i < 2; ++i)
            data.c[i] = 2 * c[i] / tau;
        vector<GridFunction1D> sol(2);
        for (int i = 0; i < 2; ++i)
            sol[i].set_grid(grid);
        vector<TimeGridFunction1D> res(2, TimeGridFunction1D(grid, tnum));

        // initial conditions
        for (int i = 0; i < 2; ++i) {
            for (int n = 0; n <= K0; ++n) {
                sol[i](0, n) = 0;
                res[i](0, 0, n) = sol[i](0, n);
            }
        }

        data.w[0][0] = data.w[0][1] = wfun(0.0);
        for (int m = 1; m <= tnum; ++m) {
            double t = m * tau, t_prev = (m - 1) * tau;
            for (int n = 0; n <= K0; ++n) {
                double x = grid.coord(0, n);
                for (int i = 0; i < 2; ++i)
                    data.g[i](0, n) = gfun[i](x, t_prev);
            }
            data.c[0] = - data.c[0];
            for (int n = 0; n <= K0; ++n) {
                double x = grid.coord(0, n);
                for (int i = 0; i < 2; ++i) {
                    if (i == 0) {
                        data.g[i](0, n) = gfun[i](x, t)
                            - OperatorValue1D(data, sol, i, 0, n);
                    }
                    else
                        data.g[i](0, n) = 0.0;
                }
            }
            data.c[0] = - data.c[0];
            data.w[0][0] = data.w[0][1] = wfun(t);
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
