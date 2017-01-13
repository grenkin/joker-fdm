/*
A simple linear conjugation problem in 1D
-u''(x) = -2 on (0, 1), -2u''(x) = -2 on (1, 3)
u(0) = 0, 2u'(3) + u(3) = 11, u'(1 - 0) = 2u'(1 + 0), u(1 - 0) = u(1 + 0)
Exact solution: u(x) = x^2 on (0, 1), u(x) = (x^2 + 1)/2 on (1, 3)
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
    return -2;
}

double exact(double x)
{
    if (x < 1)
        return x * x;
    else
        return (x * x + 1) / 2;
}

int main() {
    int K0 = 10;
    double rms_old;
    for (int test = 0; test <= 5; ++test) {
        vector<double> L(2);
        L[0] = 1;  L[1] = 2;
        vector<int> K(2);
        K[0] = K0;  K[1] = K0;
        Grid1D grid(L, K);
        Data1D data(1, grid);
        data.a[0][0] = 1;  data.a[0][1] = 2;
        data.b[0][0] = INFINITY;  data.w[0][0] = 0;
        data.b[0][1] = 1;  data.w[0][1] = 11;
        data.G[0][0] = INFINITY;
        data.f[0][0][0] = data.df[0][0][0] = zero;
        data.f[0][1][0] = data.df[0][1][0] = zero;
        vector<GridFunction1D> sol(1);
        sol[0].set_grid(grid);
        for (int j = 0; j < 2; ++j) {
            for (int n = 0; n <= K0; ++n) {
                data.g[0](j, n) = gfun(grid.coord(j, n));
                sol[0](j, n) = 0;
            }
        }
        Parameters param;
        param.max_Newton_iterations = 1;
        SolveBVP1D(data, param, sol);
        double rms = 0;
        for (int j = 0; j < 2; ++j) {
            for (int n = 0; n <= K0; ++n)
                rms += pow(sol[0](j, n) - exact(grid.coord(j, n)), 2);
        }
        rms = sqrt(rms / (2 * (K0 + 1)));
        cout << "h = " << grid.h[0] << endl;
        cout << "rms = " << rms << endl;
        if (test > 0)
            cout << "rms_old / rms = " << rms_old / rms << "\n\n";
        rms_old = rms;
        K0 *= 2;
    }

    return 0;
}
