/*
A simple linear problem in 1D
-u''(x) = -2 on (0, 1)
u(0) = 0, u'(1) + 2 * (u(1) - 2) = 0
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
    return -2;
}

double exact(double x)
{
    return x * x;
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
        data.b[0][0] = INFINITY;  data.w[0][0] = 0;
        data.b[0][1] = 2;  data.w[0][1] = 4;
        data.f[0][0][0] = data.df[0][0][0] = zero;
        vector<GridFunction1D> sol(1);
        sol[0].set_grid(grid);
        for (int n = 0; n <= K0; ++n) {
            data.g[0](0, n) = gfun(grid.coord(0, n));
            sol[0](0, n) = 0;
        }
        Parameters1D param;
        param.max_Newton_iterations = 1;
        SolveBVP1D(data, param, sol);
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
