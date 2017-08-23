/*
A linear problem in 2D
-2\Delta u = -6y, 0 < y < 1
-\Delta u = -6y, 1 < y < 2
-2u_x(0, y) = -2, 0 < y < 1
-u_x(0, y) = -1, 1 < y < 2
2u_x(1, y) = 2, 0 < y < 1
u_x(1, y) = 1, 1 < y < 2
-2u_y(x, 0) + 2u(x, 0) = 2x
u_y(x, 2) + 3u(x, 2) = 3x + 36
2u_y(x, 1 - 0) = u_y(x, 1 + 0) = 6(u(x, 1 + 0) - u(x, 1 - 0))
Exact solution: u(x, y) = x + y^3/2, 0 < y < 1
                u(x, y) = x + y^3, 1 < y < 2
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp2d.h>

using namespace std;

double gfun(double x, double y)
{
    return -6 * y;
}

double zero(double u)
{
    return 0;
}

double exact1(double x, double y)
{
    return x + pow(y, 3) / 2;
}

double exact2(double x, double y)
{
    return x + pow(y, 3);
}

int main() {
    int K0 = 10;
    double rms_old;
    for (int test = 0; test <= 5; ++test) {
        array<vector<double>, 2> L;
        array<vector<int>, 2> K;
        L[VAR_X].resize(1);  L[VAR_Y].resize(2);
        L[VAR_X][0] = 1;
        L[VAR_Y][0] = L[VAR_Y][1] = 1;
        K[VAR_X].resize(1);  K[VAR_Y].resize(2);
        K[VAR_X][0] = 2 * K0;
        K[VAR_Y][0] = K[VAR_Y][1] = 3 * K0 / 2;
        Grid2D grid(L, K);
        Data2D data(1, grid);
        data.a[0][0][0] = 2;
        data.a[0][0][1] = 1;
        for (int jY = 0; jY < 2; ++jY) {
            for (int nY = 0; nY <= grid.K[VAR_Y][jY]; ++nY) {
                data.b[0](VAR_X, 0, jY, nY) = 0;
                data.w[0](VAR_X, 0, jY, nY) = (jY == 0 ? -2 : -1);
                data.b[0](VAR_X, 1, jY, nY) = 0;
                data.w[0](VAR_X, 1, jY, nY) = (jY == 0 ? 2 : 1);
            }
        }
        for (int nX = 0; nX <= grid.K[VAR_X][0]; ++nX) {
            data.b[0](VAR_Y, 0, 0, nX) = 2;
            data.w[0](VAR_Y, 0, 0, nX) = 2 * grid.coordX(0, nX);
            data.b[0](VAR_Y, 1, 0, nX) = 3;
            data.w[0](VAR_Y, 1, 0, nX) = 3 * grid.coordX(0, nX) + 36;
        }
        data.G[0][VAR_Y][0][0] = 6;
        data.f[0][0][0][0] = data.df[0][0][0][0] = zero;
        data.f[0][0][1][0] = data.df[0][0][1][0] = zero;
        vector<GridFunction2D> sol(1);
        sol[0].set_grid(grid);
        for (int nX = 0; nX <= grid.K[VAR_X][0]; ++nX) {
            for (int jY = 0; jY < 2; ++jY) {
                for (int nY = 0; nY <= grid.K[VAR_Y][jY]; ++nY) {
                    data.g[0](0, jY, nX, nY) =
                        gfun(grid.coordX(0, nX), grid.coordY(jY, nY));
                    sol[0](0, jY, nX, nY) = 0;
                }
            }
        }
        Parameters2D param;
        param.max_Newton_iterations = 1;
        SolveBVP2D(data, param, sol);
        double rms = 0;
        for (int nX = 0; nX <= grid.K[VAR_X][0]; ++nX) {
            for (int nY = 0; nY <= grid.K[VAR_Y][0]; ++nY) {
                rms += pow(sol[0](0, 0, nX, nY) - exact1(
                    grid.coordX(0, nX), grid.coordY(0, nY)), 2);
            }
        }
        for (int nX = 0; nX <= grid.K[VAR_X][0]; ++nX) {
            for (int nY = 0; nY <= grid.K[VAR_Y][1]; ++nY) {
                rms += pow(sol[0](0, 1, nX, nY) - exact2(
                    grid.coordX(0, nX), grid.coordY(1, nY)), 2);
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
