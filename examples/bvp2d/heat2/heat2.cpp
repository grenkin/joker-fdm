#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp2d.h>

using namespace std;

double Tmax = 773;
double StBol = 5.67e-8;
double a = 0.0515;
double beta = 1;
double kappa = 10;
double alpha = 1 / (3 * kappa);
double n = 1;
double b = 4 * StBol * pow(Tmax, 3) * beta * pow(n, 2);
double c = 10;
double gamma = 0.3;

double phi_phi(double u)
{
    return beta * u;
}

double d_phi_phi(double u)
{
    return beta;
}

double phi_theta(double u)
{
    return - beta * pow(u, 4);
}

double d_phi_theta(double u)
{
    return - beta * 4 * pow(u, 3);
}

double theta_theta(double u)
{
    return b * pow(u, 4);
}

double d_theta_theta(double u)
{
    return b * 4 * pow(u, 3);
}

double theta_phi(double u)
{
    return - b * u;
}

double d_theta_phi(double u)
{
    return - b;
}

int main() {
    double Len = 1;
    int K0 = 300;
    array<vector<double>, 2> L;
    array<vector<int>, 2> K;
    L[VAR_X].resize(1);  L[VAR_Y].resize(1);
    L[VAR_X][0] = Len;
    L[VAR_Y][0] = Len;
    K[VAR_X].resize(1);  K[VAR_Y].resize(1);
    K[VAR_X][0] = K0;
    K[VAR_Y][0] = K0;
    Grid2D grid(L, K);
    Data2D data(2, grid);
    data.a[0][0][0] = alpha;
    data.a[1][0][0] = a;
    for (int nY = 0; nY <= grid.K[VAR_Y][0]; ++nY) {
        double thetab = 0.5 + 0.5 * grid.coordY(0, nY) / Len;
        data.b[0](VAR_X, 0, 0, nY) = gamma;
        data.w[0](VAR_X, 0, 0, nY) = gamma * pow(thetab, 4);
        data.b[1](VAR_X, 0, 0, nY) = c;
        data.w[1](VAR_X, 0, 0, nY) = c * thetab;

        data.b[0](VAR_X, 1, 0, nY) = gamma;
        data.w[0](VAR_X, 1, 0, nY) = gamma * pow(thetab, 4);
        data.b[1](VAR_X, 1, 0, nY) = c;
        data.w[1](VAR_X, 1, 0, nY) = c * thetab;
    }
    for (int nX = 0; nX <= grid.K[VAR_X][0]; ++nX) {
        double thetab = 0.5;
        data.b[0](VAR_Y, 0, 0, nX) = gamma;
        data.w[0](VAR_Y, 0, 0, nX) = gamma * pow(thetab, 4);
        data.b[1](VAR_Y, 0, 0, nX) = c;
        data.w[1](VAR_Y, 0, 0, nX) = c * thetab;

        thetab = 1;
        data.b[0](VAR_Y, 1, 0, nX) = gamma;
        data.w[0](VAR_Y, 1, 0, nX) = gamma * pow(thetab, 4);
        data.b[1](VAR_Y, 1, 0, nX) = c;
        data.w[1](VAR_Y, 1, 0, nX) = c * thetab;
    }
    data.f[0][0][0][0] = phi_phi;
    data.df[0][0][0][0] = d_phi_phi;
    data.f[0][0][0][1] = phi_theta;
    data.df[0][0][0][1] = d_phi_theta;
    data.f[1][0][0][0] = theta_phi;
    data.df[1][0][0][0] = d_theta_phi;
    data.f[1][0][0][1] = theta_theta;
    data.df[1][0][0][1] = d_theta_theta;
    vector<GridFunction2D> sol(2);
    sol[0].set_grid(grid);
    sol[1].set_grid(grid);
    for (int nX = 0; nX <= grid.K[VAR_X][0]; ++nX) {
        for (int nY = 0; nY <= grid.K[VAR_Y][0]; ++nY) {
            for (int i = 0; i < 2; ++i) {
                data.g[i](0, 0, nX, nY) = 0;
                sol[i](0, 0, nX, nY) = 0;
            }
        }
    }
    SolveBVP2D(data, Parameters2D(), sol);

    ofstream fout("output.txt");
    for (int nX = 0; nX <= grid.K[VAR_X][0]; ++nX) {
        for (int nY = 0; nY <= grid.K[VAR_Y][0]; ++nY) {
            double x = grid.coordX(0, nX);
            double y = grid.coordY(0, nY);
            fout << x << "  " << y << "  "
                << sol[0](0, 0, nX, nY) << "  "
                << sol[1](0, 0, nX, nY) << "\n";
        }
        fout << "\n";
    }

    return 0;
}
