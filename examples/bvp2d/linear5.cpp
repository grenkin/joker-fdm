/*
A linear problem in 2D
-\Delta u = 2 * sin(x) * sin(y)
(0, pi) \times (0, pi)
u = 0 on the boundary
Exact solution: u(x, y) = sin(x) * sin(y)
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp2d.h>

using namespace std;

double gfun(double x, double y)
{
    return 2 * sin(x) * sin(y);
}

double zero(double u)
{
    return 0;
}

double exact(double x, double y)
{
    return sin(x) * sin(y);
}

int main() {
    double L0 = 3.14159265358979323846;
    int K0 = 10;
    double rms_old;
    for (int test = 0; test <= 5; ++test) {
        array<vector<double>, 2> L;
        array<vector<int>, 2> K;
        L[VAR_X].resize(1);  L[VAR_Y].resize(1);
        L[VAR_X][0] = L[VAR_Y][0] = L0;
        K[VAR_X].resize(1);  K[VAR_Y].resize(1);
        K[VAR_X][0] = K[VAR_Y][0] = K0;
        Grid2D grid(L, K);
        Data2D data(1, grid);
        data.a[0][0][0] = 1;
        for (int nY = 0; nY <= grid.K[VAR_Y][0]; ++nY) {
            for (int s = 0; s < 2; ++s) {
                data.b[0](VAR_X, s, 0, nY) = INFINITY;
                data.w[0](VAR_X, s, 0, nY) = 0;
            }
        }
        for (int nX = 0; nX <= grid.K[VAR_X][0]; ++nX) {
            for (int s = 0; s < 2; ++s) {
                data.b[0](VAR_Y, s, 0, nX) = INFINITY;
                data.w[0](VAR_Y, s, 0, nX) = 0;
            }
        }
        data.f[0][0][0][0] = data.df[0][0][0][0] = zero;
        vector<GridFunction2D> sol(1);
        sol[0].set_grid(grid);
        for (int nX = 0; nX <= K0; ++nX) {
            for (int nY = 0; nY <= K0; ++nY) {
                data.g[0](0, 0, nX, nY) =
                    gfun(grid.coordX(0, nX), grid.coordY(0, nY));
                sol[0](0, 0, nX, nY) = 0;
            }
        }
        Parameters2D param;
        param.max_Newton_iterations = 1;
        SolveBVP2D(data, param, sol);
        double rms = 0;
        for (int nX = 0; nX <= K0; ++nX) {
            for (int nY = 0; nY <= K0; ++nY) {
                rms += pow(sol[0](0, 0, nX, nY) - exact(
                    grid.coordX(0, nX), grid.coordY(0, nY)), 2);
            }
        }
        rms = sqrt(rms / pow(K0 + 1, 2));
        cout << "h = " << grid.h[VAR_X][0] << endl;
        cout << "rms = " << rms << endl;
        if (test > 0)
            cout << "rms_old / rms = " << rms_old / rms << "\n\n";
        rms_old = rms;
        K0 *= 2;
    }
    return 0;
}
