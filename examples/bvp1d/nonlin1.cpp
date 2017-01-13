/*
A nonlinear problem in 1D
-u''(x) + u^2 = -2 + x^4 on (0, 1)
u'(0) = 0, u(1) = 1
Exact solution: u(x) = x ^ 2
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp1d.h>

using namespace std;

double zero(double x)
{
    return 0;
}

double gfun(double x)
{
    return -2 + pow(x, 4);
}

double exact(double x)
{
    return x * x;
}

double f(double u)
{
    return u * u;
}

double df(double u)
{
    return 2 * u;
}

int main() {
    double L0 = 1;
    int K0 = 10;
    double rms_old;
    for (int test = 0; test <= 5; ++test) {
        vector<double> L(1);  L[0] = L0;
        vector<int> K(1);  K[0] = K0;
        Grid1D grid(L, K);
        Data1D data(1, grid);
        data.a[0][0] = 1;
        data.b[0][0] = 0;  data.w[0][0] = 0;
        data.b[0][1] = INFINITY;  data.w[0][1] = 1;
        data.f[0][0][0] = f;  data.df[0][0][0] = df;
        vector<GridFunction1D> sol(1);
        sol[0].set_grid(grid);
        for (int n = 0; n <= K0; ++n) {
            data.g[0](0, n) = gfun(grid.coord(0, n));
            sol[0](0, n) = 0;
        }
        SolveBVP1D(data, Parameters(), sol);
        double rms = 0;
        for (int n = 0; n <= K0; ++n)
            rms += pow(sol[0](0, n) - exact(grid.coord(0, n)), 2);
        rms = sqrt(rms / (K0 + 1));
        cout << "h = " << grid.h[0] << endl;
        cout << "rms = " << rms << endl;
        if (test > 0)
            cout << "rms_old / rms = " << rms_old / rms << "\n\n";
        rms_old = rms;
        K0 *= 2;
    }

    return 0;
}
