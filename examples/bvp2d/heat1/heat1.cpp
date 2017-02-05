#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <joker-fdm/bvp2d.h>

using namespace std;

double Tmax = 773;
double StBol = 5.67e-8;
// 1 - glass, 2 - air
double a1 = 1.2;
double beta1 = 100;
double kappa1 = 100;
double alpha1 = 1 / (3 * kappa1);
double n1 = 1.47;
double b1 = 4 * StBol * pow(Tmax, 3) * beta1 * pow(n1, 2);
double a2 = 0.0515;
double beta2 = 1;
double kappa2 = 10;
double alpha2 = 1 / (3 * kappa2);
double n2 = 1;
double b2 = 4 * StBol * pow(Tmax, 3) * beta2 * pow(n2, 2);
double c = 10;
double gamma = 0.3;

double phi_phi_1(double u)
{
    return pow(n1, 2) * beta1 * u;
}

double d_phi_phi_1(double u)
{
    return pow(n1, 2) * beta1;
}

double phi_phi_2(double u)
{
    return pow(n2, 2) * beta2 * u;
}

double d_phi_phi_2(double u)
{
    return pow(n2, 2) * beta2;
}

double phi_theta_1(double u)
{
    return - pow(n1, 2) * beta1 * pow(u, 4);
}

double d_phi_theta_1(double u)
{
    return - pow(n1, 2) * beta1 * 4 * pow(u, 3);
}

double phi_theta_2(double u)
{
    return - pow(n2, 2) * beta2 * pow(u, 4);
}

double d_phi_theta_2(double u)
{
    return - pow(n2, 2) * beta2 * 4 * pow(u, 3);
}

double theta_theta_1(double u)
{
    return b1 * pow(u, 4);
}

double d_theta_theta_1(double u)
{
    return b1 * 4 * pow(u, 3);
}

double theta_theta_2(double u)
{
    return b2 * pow(u, 4);
}

double d_theta_theta_2(double u)
{
    return b2 * 4 * pow(u, 3);
}

double theta_phi_1(double u)
{
    return - b1 * u;
}

double d_theta_phi_1(double u)
{
    return - b1;
}

double theta_phi_2(double u)
{
    return - b2 * u;
}

double d_theta_phi_2(double u)
{
    return - b2;
}

struct n_params {
    double n1, n2;
};

double psi_ij (double nu, double n_ij)
{
    double radicand = 1 - pow(n_ij, 2) * (1 - pow(nu, 2));
    if (radicand > 0)
        return sqrt(radicand);
    else
        return 0;
}

double R_ij (double nu, double n_ij)
{
    if (nu == 0)
        return 1;
    else {
        double frac1 = (psi_ij(nu, n_ij) - n_ij * nu)
            / (psi_ij(nu, n_ij) + n_ij * nu);
        double frac2 = (n_ij * psi_ij(nu, n_ij) - nu)
            / (n_ij * psi_ij(nu, n_ij) + nu);
        return 0.5 * (pow(frac1, 2) + pow(frac2, 2));
    }
}

double integrand1 (double nu, void* params)
{
    n_params* par = (n_params*)params;
    return nu * (1 - R_ij(nu, par->n1 / par->n2));
}

double integrand2 (double nu, void* params)
{
    n_params* par = (n_params*)params;
    return pow(nu, 2) * (R_ij(nu, par->n1 / par->n2) + R_ij(nu, par->n2 / par->n1));
}

double calc_G (double n1, double n2)
{
    double int1, int2, error;
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &integrand1;
    n_params par;  par.n1 = n1;  par.n2 = n2;
    F.params = &par;
    gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000, w, &int1, &error);
    F.function = &integrand2;
    gsl_integration_qags(&F, 0, 1, 0, 1e-7, 1000, w, &int2, &error);
    gsl_integration_workspace_free(w);
    return pow(n1, 2) * int1 / (3 * int2);
}

int main() {
    double Len = 0.25, d = 0.1, d1 = 0.1, d2 = 0.1;
    int K1 = 200, K2 = 200, K3 = 100;
    array<vector<double>, 2> L;
    array<vector<int>, 2> K;
    L[VAR_X].resize(3);  L[VAR_Y].resize(3);
    L[VAR_X][0] = d1;  L[VAR_X][1] = d;  L[VAR_X][2] = Len - (d1 + d);
    L[VAR_Y][0] = d2;  L[VAR_Y][1] = d;  L[VAR_Y][2] = Len - (d2 + d);
    K[VAR_X].resize(3);  K[VAR_Y].resize(3);
    K[VAR_X][0] = K1;  K[VAR_X][1] = K2;  K[VAR_X][2] = K3;
    K[VAR_Y][0] = K1;  K[VAR_Y][1] = K2;  K[VAR_Y][2] = K3;
    Grid2D grid(L, K);
    Data2D data(2, grid);
    for (int jX = 0; jX < 3; ++jX) {
        for (int jY = 0; jY < 3; ++jY) {
            if (jX == 1 && jY == 1) {
                data.a[0][jX][jY] = pow(n2, 2) * alpha2;
                data.a[1][jX][jY] = a2;
            }
            else {
                data.a[0][jX][jY] = pow(n1, 2) * alpha1;
                data.a[1][jX][jY] = a1;
            }
        }
    }
    for (int jY = 0; jY < 3; ++jY) {
        for (int nY = 0; nY <= grid.K[VAR_Y][jY]; ++nY) {
            double thetab = 0.5 + 0.5 * grid.coordY(jY, nY) / Len;
            data.b[0](VAR_X, 0, jY, nY) = gamma * pow(n1, 2);
            data.w[0](VAR_X, 0, jY, nY) = gamma * pow(n1, 2) * pow(thetab, 4);
            data.b[1](VAR_X, 0, jY, nY) = c;
            data.w[1](VAR_X, 0, jY, nY) = c * thetab;

            data.b[0](VAR_X, 1, jY, nY) = gamma * pow(n1, 2);
            data.w[0](VAR_X, 1, jY, nY) = gamma * pow(n1, 2) * pow(thetab, 4);
            data.b[1](VAR_X, 1, jY, nY) = c;
            data.w[1](VAR_X, 1, jY, nY) = c * thetab;
        }
    }
    for (int jX = 0; jX < 3; ++jX) {
        for (int nX = 0; nX <= grid.K[VAR_X][jX]; ++nX) {
            double thetab = 0.5;
            data.b[0](VAR_Y, 0, jX, nX) = gamma * pow(n1, 2);
            data.w[0](VAR_Y, 0, jX, nX) = gamma * pow(n1, 2) * pow(thetab, 4);
            data.b[1](VAR_Y, 0, jX, nX) = c;
            data.w[1](VAR_Y, 0, jX, nX) = c * thetab;

            thetab = 1;
            data.b[0](VAR_Y, 1, jX, nX) = gamma * pow(n1, 2);
            data.w[0](VAR_Y, 1, jX, nX) = gamma * pow(n1, 2) * pow(thetab, 4);
            data.b[1](VAR_Y, 1, jX, nX) = c;
            data.w[1](VAR_Y, 1, jX, nX) = c * thetab;
        }
    }
    double INF = 1e6;
    for (int i = 0; i < 2; ++i) {
        for (int jX = 0; jX < 2; ++jX) {
            for (int jY = 0; jY < 3; ++jY)
                data.G[i][VAR_X][jX][jY] = INF;
        }
        for (int jX = 0; jX < 3; ++jX) {
            for (int jY = 0; jY < 2; ++jY)
                data.G[i][VAR_Y][jX][jY] = INF;
        }
    }
    for (int jX = 0; jX < 3; ++jX) {
        for (int jY = 0; jY < 3; ++jY) {
            if (jX == 1 && jY == 1) {
                data.f[0][jX][jY][0] = phi_phi_2;
                data.df[0][jX][jY][0] = d_phi_phi_2;
                data.f[0][jX][jY][1] = phi_theta_2;
                data.df[0][jX][jY][1] = d_phi_theta_2;
                data.f[1][jX][jY][0] = theta_phi_2;
                data.df[1][jX][jY][0] = d_theta_phi_2;
                data.f[1][jX][jY][1] = theta_theta_2;
                data.df[1][jX][jY][1] = d_theta_theta_2;
            }
            else {
                data.f[0][jX][jY][0] = phi_phi_1;
                data.df[0][jX][jY][0] = d_phi_phi_1;
                data.f[0][jX][jY][1] = phi_theta_1;
                data.df[0][jX][jY][1] = d_phi_theta_1;
                data.f[1][jX][jY][0] = theta_phi_1;
                data.df[1][jX][jY][0] = d_theta_phi_1;
                data.f[1][jX][jY][1] = theta_theta_1;
                data.df[1][jX][jY][1] = d_theta_theta_1;
            }
        }
    }
    vector<GridFunction2D> sol(2);
    sol[0].set_grid(grid);
    sol[1].set_grid(grid);
    for (int jX = 0; jX < 3; ++jX) {
        for (int nX = 0; nX <= grid.K[VAR_X][jX]; ++nX) {
            for (int jY = 0; jY < 3; ++jY) {
                for (int nY = 0; nY <= grid.K[VAR_Y][jY]; ++nY) {
                    for (int i = 0; i < 2; ++i) {
                        data.g[i](jX, jY, nX, nY) = 0;
                        sol[i](jX, jY, nX, nY) = 0;
                    }
                }
            }
        }
    }
    cout << "Calculate without reflection and refraction...\n";
    SolveBVP2D(data, Parameters2D(), sol);
    vector<GridFunction2D> sol_wo = sol;
    double G0 = calc_G(n1, n2);
    data.G[0][VAR_X][0][1] = data.G[0][VAR_X][1][1] = G0;
    data.G[0][VAR_Y][1][0] = data.G[0][VAR_Y][1][1] = G0;
    cout << "Calculate with reflection and refraction...\n";
    SolveBVP2D(data, Parameters2D(), sol);

    ofstream fout1("output_without.txt");
    ofstream fout2("output_with.txt");
    for (int jX = 0; jX < 3; ++jX) {
        for (int nX = 0; nX <= grid.K[VAR_X][jX]; ++nX) {
            for (int jY = 0; jY < 3; ++jY) {
                for (int nY = 0; nY <= grid.K[VAR_Y][jY]; ++nY) {
                    double x = grid.coordX(jX, nX);
                    double y = grid.coordY(jY, nY);
                    if (nX == 0 && jX > 0)
                        x += grid.h[VAR_X][jX] / 100;
                    if (nY == 0 && jY > 0)
                        y += grid.h[VAR_Y][jY] / 100;
                    fout1 << x << "  " << y << "  "
                        << sol_wo[0](jX, jY, nX, nY) << "  "
                        << sol_wo[1](jX, jY, nX, nY) << "\n";
                    fout2 << x << "  " << y << "  "
                        << sol[0](jX, jY, nX, nY) << "  "
                        << sol[1](jX, jY, nX, nY) << "\n";
                }
            }
            fout1 << "\n";
            fout2 << "\n";
        }
    }
    return 0;
}
