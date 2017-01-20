/*
A linear problem in 2D
-\Delta u + v = 2 * sin(x) * sin(y) + sin(x)
-\Delta v = sin(x)
-u_x + u|_{x=0} = -sin(y)
u_x + u|_{x=pi} = -sin(y)
-u_y + u|_{y=0} = -sin(x)
u_y + u|_{y=pi} = -sin(x)
-v_x + v|_{x=0} = -1
v_x + v|_{x=pi} = -1
-v_y + v|_{y=0} = sin(x)
v_y + v|_{y=pi} = sin(x)
Exact solution: u(x, y) = sin(x) * sin(y), v(x, y) = sin(x)
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp2d.h>

using namespace std;

double gfun1(double x, double y)
{
    return 2 * sin(x) * sin(y) + sin(x);
}

double gfun2(double x, double y)
{
    return sin(x);
}

double zero(double u)
{
    return 0;
}

double ident(double u)
{
    return u;
}

double one(double u)
{
    return 1;
}

double exact1(double x, double y)
{
    return sin(x) * sin(y);
}

double exact2(double x, double y)
{
    return sin(x);
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
        Data2D data(2, grid);
        data.a[0][0][0] = data.a[1][0][0] = 1;
        for (int nY = 0; nY <= grid.K[VAR_Y][0]; ++nY) {
            for (int s = 0; s < 2; ++s) {
                data.b[0](VAR_X, s, 0, nY) = 1;
                data.w[0](VAR_X, s, 0, nY) = -sin(grid.coordY(0, nY));
                data.b[1](VAR_X, s, 0, nY) = 1;
                data.w[1](VAR_X, s, 0, nY) = -1;
            }
        }
        for (int nX = 0; nX <= grid.K[VAR_X][0]; ++nX) {
            for (int s = 0; s < 2; ++s) {
                data.b[0](VAR_Y, s, 0, nX) = 1;
                data.w[0](VAR_Y, s, 0, nX) = -sin(grid.coordX(0, nX));
                data.b[1](VAR_Y, s, 0, nX) = 1;
                data.w[1](VAR_Y, s, 0, nX) = sin(grid.coordX(0, nX));
            }
        }
        for (int i = 0; i < 2; ++i) {
            for (int k = 0; k < 2; ++k)
                data.f[i][0][0][k] = data.df[i][0][0][k] = zero;
        }
        data.f[0][0][0][1] = ident;  data.df[0][0][0][1] = one;
        vector<GridFunction2D> sol(2);
        sol[0].set_grid(grid);
        sol[1].set_grid(grid);
        for (int nX = 0; nX <= K0; ++nX) {
            for (int nY = 0; nY <= K0; ++nY) {
                data.g[0](0, 0, nX, nY) =
                    gfun1(grid.coordX(0, nX), grid.coordY(0, nY));
                data.g[1](0, 0, nX, nY) =
                    gfun2(grid.coordX(0, nX), grid.coordY(0, nY));
                sol[0](0, 0, nX, nY) = 0;
                sol[1](0, 0, nX, nY) = 0;
            }
        }
        Parameters2D param;
        param.max_Newton_iterations = 1;
        SolveBVP2D(data, param, sol);
        double rms = 0;
        for (int nX = 0; nX <= K0; ++nX) {
            for (int nY = 0; nY <= K0; ++nY) {
                rms += pow(sol[0](0, 0, nX, nY) - exact1(
                    grid.coordX(0, nX), grid.coordY(0, nY)), 2)
                    + pow(sol[1](0, 0, nX, nY) - exact2(
                    grid.coordX(0, nX), grid.coordY(0, nY)), 2);
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
