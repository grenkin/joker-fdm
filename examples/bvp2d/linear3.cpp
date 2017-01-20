/*
A linear problem in 2D
-\Delta u = -2, 0 < x < 1
-2\Delta u = -2, 1 < x < 2
-u_x(0, y) + 5u(0, y) = 5y
2u_x(2, y) + 4u(2, y) = 4y + 16
-u_y(x, 0) = -1, 0 < x < 1
-2u_y(x, 0) = -2, 1 < x < 2
u_y(x, 1) = 1, 0 < x < 1
2u_y(x, 1) = 2, 1 < x < 2
u_x(1 - 0, y) = 2u_x(1 + 0, y) = 4(u(1 + 0, y) - u(1 - 0, y))
Exact solution: u(x, y) = x^2 + y, 0 < x < 1
                u(x, y) = x^2/2 + y + 1, 1 < x < 2
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp2d.h>

using namespace std;

double gfun(double x, double y)
{
    return -2;
}

double zero(double u)
{
    return 0;
}

double exact1(double x, double y)
{
    return x * x + y;
}

double exact2(double x, double y)
{
    return x * x / 2 + y + 1;
}

int main() {
    int K0 = 10;
    double rms_old;
    for (int test = 0; test <= 5; ++test) {
        array<vector<double>, 2> L;
        array<vector<int>, 2> K;
        L[VAR_X].resize(2);  L[VAR_Y].resize(1);
        L[VAR_X][0] = L[VAR_X][1] = 1;
        L[VAR_Y][0] = 1;
        K[VAR_X].resize(2);  K[VAR_Y].resize(1);
        K[VAR_X][0] = K[VAR_X][1] = K0;
        K[VAR_Y][0] = K0;
        Grid2D grid(L, K);
        Data2D data(1, grid);
        data.a[0][0][0] = 1;
        data.a[0][1][0] = 2;
        for (int nY = 0; nY <= grid.K[VAR_Y][0]; ++nY) {
            data.b[0](VAR_X, 0, 0, nY) = 5;
            data.w[0](VAR_X, 0, 0, nY) = 5 * grid.coordY(0, nY);
            data.b[0](VAR_X, 1, 0, nY) = 4;
            data.w[0](VAR_X, 1, 0, nY) = 4 * grid.coordY(0, nY) + 16;
        }
        for (int jX = 0; jX < 2; ++jX) {
            for (int nX = 0; nX <= grid.K[VAR_X][jX]; ++nX) {
                data.b[0](VAR_Y, 0, jX, nX) = 0;
                data.w[0](VAR_Y, 0, jX, nX) = (jX == 0 ? -1 : -2);
                data.b[0](VAR_Y, 1, jX, nX) = 0;
                data.w[0](VAR_Y, 1, jX, nX) = (jX == 0 ? 1 : 2);
            }
        }
        data.G[0][VAR_X][0][0] = 4;
        data.f[0][0][0][0] = data.df[0][0][0][0] = zero;
        data.f[0][1][0][0] = data.df[0][1][0][0] = zero;
        vector<GridFunction2D> sol(1);
        sol[0].set_grid(grid);
        for (int jX = 0; jX < 2; ++jX) {
            for (int nX = 0; nX <= K0; ++nX) {
                for (int nY = 0; nY <= K0; ++nY) {
                    data.g[0](jX, 0, nX, nY) =
                        gfun(grid.coordX(jX, nX), grid.coordY(0, nY));
                    sol[0](jX, 0, nX, nY) = 0;
                }
            }
        }
        Parameters2D param;
        param.max_Newton_iterations = 1;
        SolveBVP2D(data, param, sol);
        double rms = 0;
        for (int nX = 0; nX <= K0; ++nX) {
            for (int nY = 0; nY <= K0; ++nY) {
                rms += pow(sol[0](0, 0, nX, nY) - exact1(
                    grid.coordX(0, nX), grid.coordY(0, nY)), 2);
            }
        }
        for (int nX = 0; nX <= K0; ++nX) {
            for (int nY = 0; nY <= K0; ++nY) {
                rms += pow(sol[0](1, 0, nX, nY) - exact2(
                    grid.coordX(1, nX), grid.coordY(0, nY)), 2);
            }
        }
        rms = sqrt(rms / (2 * pow(K0 + 1, 2)));
        cout << "h = " << grid.h[VAR_X][0] << endl;
        cout << "rms = " << rms << endl;
        if (test > 0)
            cout << "rms_old / rms = " << rms_old / rms << "\n\n";
        rms_old = rms;
        K0 *= 2;
    }
    return 0;
}
