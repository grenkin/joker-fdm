/*
An initial-boundary value problem in 1D
0.5u_t - 0.5u_xx = 2*pi^2*sin(2*pi*(x-t)) - pi*cos(2*pi*(x-t)), x \in (0, 1)
u(0, t) = -sin(2*pi*t),  u(1, t) = sin(2*pi*(1-t))
u(x, 0) = sin(2*pi*x)
Exact solution: u(x, t) = sin(2*pi*(x-t))
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

double exact (double x, double t)
{
    return sin(2 * pi * (x - t));
}

double gfun (double x, double t)
{
    return 2 * pi * pi * sin(2 * pi * (x - t)) - pi * cos(2 * pi * (x - t));
}

double w1fun (double t)
{
    return - sin(2 * pi * t);
}

double w2fun (double t)
{
    return sin(2 * pi * (1 - t));
}

double initfun (double x)
{
    return sin(2 * pi * x);
}

int main() {
    double L0 = 1;
    double T = 30;
    int K0 = 10;
    int tnum = 10;
    double rms_old;
    for (int test = 0; test <= 5; ++test) {
        vector<double> L(1);  L[0] = L0;
        vector<int> K(1);  K[0] = K0;
        Grid1D grid(L, K);
        Data1D data(1, grid);
        data.a[0][0] = 0.5;
        data.b[0][0] = INFINITY;  // data.w[0][0] will be defined later
        data.b[0][1] = INFINITY;  // data.w[0][1] will be defined later
        data.f[0][0][0] = data.df[0][0][0] = zero;
        double c = 0.5;
        double tau = T / tnum;
        data.c[0] = c / tau;
        vector<GridFunction1D> sol(1);
        sol[0].set_grid(grid);
        TimeGridFunction1D res(grid, tnum);

        // initial condition
        for (int n = 0; n <= K0; ++n) {
            sol[0](0, n) = initfun(grid.coord(0, n));
            res(0, 0, n) = sol[0](0, n);
        }

        for (int m = 1; m <= tnum; ++m) {
            double t = m * tau;
            data.w[0][0] = w1fun(t);
            data.w[0][1] = w2fun(t);
            for (int n = 0; n <= K0; ++n) {
                data.g[0](0, n) = gfun(grid.coord(0, n), t) + c / tau * sol[0](0, n);
            }
            Parameters1D param;
            param.max_Newton_iterations = 1;
            SolveBVP1D(data, param, sol);
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
        tnum *= 4;
    }

    return 0;
}
