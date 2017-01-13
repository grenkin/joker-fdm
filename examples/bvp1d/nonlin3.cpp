/*
A nonlinear conjugation problem in 1D
-u''(x) + u^4 = sin(x) + sin^4(x) on (0, pi/2)
-4u''(x) + u = 2 * sin(x/2 + pi/4) on (pi/2, 3*pi/2)
-u'(0) + u(0) = -1, u(3*pi/2) = 0
u'(pi/2 - 0) = 4u'(pi/2 + 0)
Exact solution: u(x) = sin(x), x in (0, pi/2);
                u(x) = sin(x/2 + pi/4), x in (pi/2, 3*pi/2)
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp1d.h>

using namespace std;

const double pi = 3.14159265358979323846;

double gfun(double x)
{
    if (x < pi / 2)
        return sin(x) + pow(sin(x), 4);
    else
        return 2 * sin(x / 2 + pi / 4);
}

double exact(double x)
{
    return x < pi / 2  ?  sin(x)  :  sin(x / 2 + pi / 4);
}

double f(double u)
{
    return pow(u, 4);
}

double df(double u)
{
    return 4 * pow(u, 3);
}

double ident(double u)
{
    return u;
}

double one(double u)
{
    return 1;
}

int main() {
    int K0 = 10;
    double rms_old;
    for (int test = 0; test <= 5; ++test) {
        vector<double> L(2);  L[0] = pi / 2;  L[1] = pi;
        vector<int> K(2);  K[0] = K[1] = K0;
        Grid1D grid(L, K);
        Data1D data(1, grid);
        data.a[0][0] = 1;  data.a[0][1] = 4;
        data.b[0][0] = 1;  data.w[0][0] = -1;
        data.b[0][1] = INFINITY;  data.w[0][1] = 0;
        data.G[0][0] = INFINITY;
        data.f[0][0][0] = f;  data.df[0][0][0] = df;
        data.f[0][1][0] = ident;  data.df[0][1][0] = one;
        vector<GridFunction1D> sol(1);
        sol[0].set_grid(grid);
        for (int j = 0; j < 2; ++j) {
            for (int n = 0; n <= K0; ++n) {
                data.g[0](j, n) = gfun(grid.coord(j, n));
                sol[0](j, n) = 0;
            }
        }
        SolveBVP1D(data, Parameters(), sol);
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
