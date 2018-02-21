// A nonstationary radiative heat transfer model

#include <iostream>
#include <fstream>
#include <cmath>
#include <joker-fdm/bvp1d.h>

using namespace std;

const double pi = 3.14159265358979323846;

const double L = 0.5;  // length of the layer
const double T = 6e-8;  // length of the time interval

// parameters of air at 400 C and normal atmospheric pressure
const double k = 0.0521;  // thermal conductivity
const double rho = 0.524;  // density
const double c_p = 1068;  // heat capacity
const double n = 1.0;  // refractive index

// parameters of the boundary
const double h = 10;  // heat transfer coefficient
const double eps = 0.8;  // emissivity of the boundary

const double Tmax = 773;  // maximum temperature in the unnormalized model
const double sigma = 5.67e-8;  // Stefan-Boltzmann constant
const double a = k / (rho * c_p);
const double b = 4 * sigma * pow(n, 2) * pow(Tmax, 3) / (rho * c_p);

const double kappa = 10;  // extinction coefficient
const double kappa_a = 1;  // absorption coefficient
const double alpha = 1 / (3 * kappa);

const double beta = h / (rho * c_p);
const double gamma = eps / (2 * (2 - eps));

const double c = 3e8;  // the speed of light

const double theta_b = 0.4;

const int N = 1000;  // the number of spatial grid points
const int M = 10000;  // the number of temporal grid points

// coefficients in the initial conditions
const double A = 0.3;
const double B = 0.3;

//initial functions

double zeta_init (double x)
{
    return A * (1 - cos(2 * pi * x / L));
}

double xi_init (double x)
{
    return B * (1 - cos(2 * pi * x / L));
}

// reaction terms functions

double theta_theta(double u)
{
    return b * kappa_a * pow(u, 4);
}

double d_theta_theta(double u)
{
    return b * kappa_a * 4 * pow(u, 3);
}

double theta_phi(double u)
{
    return - b * kappa_a * u;
}

double d_theta_phi(double u)
{
    return - b * kappa_a;
}

double phi_theta(double u)
{
    return - kappa_a * pow(u, 4);
}

double d_phi_theta(double u)
{
    return - kappa_a * 4 * pow(u, 3);
}

double phi_phi(double u)
{
    return kappa_a * u;
}

double d_phi_phi(double u)
{
    return kappa_a;
}

int main() {
    vector<double> Lvec(1);  Lvec[0] = L;
    vector<int> Kvec(1);  Kvec[0] = N;
    Grid1D grid(Lvec, Kvec);
    Data1D data(2, grid);
    const double phi_b = pow(theta_b, 4);
    data.a[0][0] = a;
    data.a[1][0] = alpha;
    data.b[0][0] = data.b[0][1] = beta;
    data.w[0][0] = data.w[0][1] = beta * theta_b;
    data.b[1][0] = data.b[1][1] = gamma;
    data.w[1][0] = data.w[1][1] = gamma * phi_b;
    data.f[0][0][0] = theta_theta;  data.df[0][0][0] = d_theta_theta;
    data.f[0][0][1] = theta_phi;  data.df[0][0][1] = d_theta_phi;
    data.f[1][0][0] = phi_theta;  data.df[1][0][0] = d_phi_theta;
    data.f[1][0][1] = phi_phi;  data.df[1][0][1] = d_phi_phi;
    vector<double> cvec(2);
    cvec[0] = 1; cvec[1] = 1 / c;
    double tau = T / M;
    for (int i = 0; i < 2; ++i)
        data.c[i] = cvec[i] / tau;

    vector<GridFunction1D> sol(2);
    vector<TimeGridFunction1D> sol_time(2);
    for (int i = 0; i < 2; ++i) {
        sol[i].set_grid(grid);
        sol_time[i].set_grid(grid, M);
    }

    // initial conditions
    for (int n = 0; n <= N; ++n) {
        sol[0](0, n) = theta_b + zeta_init(grid.coord(0, n));
        sol[1](0, n) = pow(theta_b, 4) + xi_init(grid.coord(0, n));
        for (int i = 0; i < 2; ++i)
            sol_time[i](0, 0, n) = sol[i](0, n);
    }

    // solve the nonstationary problem
    for (int m = 1; m <= M; ++m) {
        cout << "m = " << m << " of " << M << endl;
        for (int n = 0; n <= N; ++n) {
            data.g[0](0, n) = cvec[0] / tau * sol[0](0, n);
            data.g[1](0, n) = cvec[1] / tau * sol[1](0, n);
        }
        Parameters1D param;
        param.max_linear_sys_iterations = 100000;
        param.linear_sys_tol = 1e-17;
        SolveBVP1D(data, param, sol);
        for (int i = 0; i < 2; ++i) {
            for (int n = 0; n <= N; ++n)
                sol_time[i](m, 0, n) = sol[i](0, n);
        }
    }

    // solve stationary problems
    Data1D data_phi(1, grid);
    data_phi.a[0][0] = alpha;
    data_phi.b[0][0] = data_phi.b[0][1] = gamma;
    data_phi.w[0][0] = data_phi.w[0][1] = gamma * phi_b;
    data_phi.f[0][0][0] = phi_phi;  data_phi.df[0][0][0] = d_phi_phi;

    vector<GridFunction1D> sol_phi(1);
    TimeGridFunction1D sol_time_phi;
    sol_phi[0].set_grid(grid);
    sol_time_phi.set_grid(grid, M);
    // initial guess
    for (int n = 0; n <= N; ++n)
        sol_phi[0](0, n) = pow(theta_b, 4);

    for (int m = 0; m <= M; ++m) {
        cout << "stat. phi: m = " << m << " of " << M << endl;
        for (int n = 0; n <= N; ++n)
            data_phi.g[0](0, n) = kappa_a * pow(sol_time[0](m, 0, n), 4);
        Parameters1D param_phi;
        param_phi.max_linear_sys_iterations = 100000;
        param_phi.linear_sys_tol = 1e-17;
        SolveBVP1D(data_phi, param_phi, sol_phi);
        for (int n = 0; n <= N; ++n)
            sol_time_phi(m, 0, n) = sol_phi[0](0, n);
    }

    // calculate the norm of the difference
    vector<double> norm_sq(M + 1), r(M + 1);
    for (int m = 0; m <= M; ++m) {
        // use the trapezoid formula
        double integral =
            0.5 * pow(sol_time[1](m, 0, 0) - sol_time_phi(m, 0, 0), 2)
            + 0.5 * pow(sol_time[1](m, 0, N) - sol_time_phi(m, 0, N), 2);
        for (int n = 1; n <= N - 1; ++n)
            integral += pow(sol_time[1](m, 0, n) - sol_time_phi(m, 0, n), 2);
        integral *= grid.h[0];
        norm_sq[m] = integral;
        r[m] = sqrt(integral / L);
    }

    // output the result
    ofstream fmid("output_mid.txt");
    fmid.precision(20);
    for (int m = 0; m <= M; ++m) {
        fmid << m * tau << "  " << sol_time[0](m, 0, N / 2) << "  "
            << sol_time[1](m, 0, N / 2) << "    "
            << sol_time_phi(m, 0, N / 2) << endl;
    }
    ofstream fdiff("output_diff.txt");
    fdiff.precision(20);
    for (int m = 0; m <= M; ++m) {
        fdiff << m * tau << "  " << norm_sq[m] << "  " << r[m] << endl;
    }

    return 0;
}
