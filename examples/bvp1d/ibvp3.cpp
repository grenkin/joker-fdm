/*
An initial-boundary value problem in 1D
u_t - u_xx = 2(x^2 - x)(2t - 1) - 4(t^2 - t), x \in (0, 1)
-u_x(0, t) + 3u(0, t) = 2(t^2 - t)
u_x(1, t) + 5u(1, t) = 2(t^2 - t)
u(x, 0) = 0
Exact solution: u(x, t) = 2(t^2 - t)(x^2 - x)
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

double exact (double x, double t)
{
    return 2 * (t * t - t) * (x * x - x);
}

double gfun (double x, double t)
{
    return 2 * (x * x - x) * (2 * t - 1) - 4 * (t * t - t);
}

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
        Data1D data(1, grid);
        data.a[0][0] = 1;
        data.b[0][0] = 3;  // data.w[0][0] will be defined later
        data.b[0][1] = 5;  // data.w[0][1] will be defined later
        data.f[0][0][0] = data.df[0][0][0] = zero;
        double c = 1;
        double tau = T / tnum;
        data.c[0] = c / tau;
        vector<GridFunction1D> sol(1);
        sol[0].set_grid(grid);
        TimeGridFunction1D res(grid, tnum);

        // initial condition
        for (int n = 0; n <= K0; ++n) {
            sol[0](0, n) = 0;
            res(0, 0, n) = sol[0](0, n);
        }

        for (int m = 1; m <= tnum; ++m) {
            double t = m * tau;
            data.w[0][0] = data.w[0][1] = 2 * (t * t - t);
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
