/*
A nonlinear initial-boundary value problem in 1D
u_t - 2u_xx + u^2 + u = 4exp(-t) + x^2*(1-x)^2*exp(-2t), x \in (0, 1)
u(0, t) = 0,  2u_x(1, t) + u(1, t) = -2exp(-t)
u(x, 0) = x*(1-x)
Exact solution: u(x, t) = x*(1-x)*exp(-t)

Crank-Nicolson scheme is used:
2 u^m / tau + Au^m = g^{m-1} + g^m + 2 u^{m-1} / tau - Au^{m-1}
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp1d.h>

using namespace std;

double f1 (double x)
{
    return x * x + x;
}

double f2 (double x)
{
    return 2 * x + 1;
}

double exact (double x, double t)
{
    return x * (1 - x) * exp(-t);
}

double gfun (double x, double t)
{
    return 4 * exp(-t) + x * x * (1 - x) * (1 - x) * exp(-2 * t);
}

double wfun (double t)
{
    return -2 * exp(-t);
}

double initfun (double x)
{
    return x * (1 - x);
}

int main ()
{
    double L0 = 1;
    double T = 15;
    int K0 = 10;
    int tnum = 10;
    double rms_old;
    for (int test = 0; test <= 5; ++test) {
        vector<double> L(1);  L[0] = L0;
        vector<int> K(1);  K[0] = K0;
        Grid1D grid(L, K);
        Data1D data(1, grid);
        data.a[0][0] = 2;
        data.b[0][0] = INFINITY;  data.w[0][0] = 0;
        data.b[0][1] = 1;  // data.w[0][1] will be defined later
        data.f[0][0][0] = f1;
        data.df[0][0][0] = f2;
        double c = 1;
        double tau = T / tnum;
        data.c[0] = 2 * c / tau;
        vector<GridFunction1D> sol(1);
        sol[0].set_grid(grid);
        TimeGridFunction1D res(grid, tnum);

        // initial condition
        for (int n = 0; n <= K0; ++n) {
            sol[0](0, n) = initfun(grid.coord(0, n));
            res(0, 0, n) = sol[0](0, n);
        }

        data.w[0][1] = wfun(0.0);
        for (int m = 1; m <= tnum; ++m) {
            double t = m * tau, t_prev = (m - 1) * tau;
            for (int n = 0; n <= K0; ++n) {
                double x = grid.coord(0, n);
                data.g[0](0, n) = gfun(x, t_prev);
            }
            data.c[0] = - data.c[0];
            for (int n = 0; n <= K0; ++n) {
                double x = grid.coord(0, n);
                data.g[0](0, n) = gfun(x, t)
                    - OperatorValue1D(data, sol, 0, 0, n);
            }
            data.c[0] = - data.c[0];
            data.w[0][1] = wfun(t);
            SolveBVP1D(data, Parameters1D(), sol);
            for (int n = 0; n <= K0; ++n)
                res(m, 0, n) = sol[0](0, n);
        }

        double rms = 0;
        for (int m = 0; m <= tnum; ++m) {
            for (int n = 0; n <= K0; ++n)
                rms += pow(res(m, 0, n) - exact(grid.coord(0, n), m * tau), 2);
        }
        rms = sqrt(rms / (K0 + 1) / (tnum + 1));

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
