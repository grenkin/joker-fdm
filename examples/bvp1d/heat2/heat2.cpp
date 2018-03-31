// An optimal control problem for a nonstationary radiative heat transfer model
//   The control parameter is the coefficient in boundary conditions for phi
//   Problem 1: minimize theta and phi
//   Problem 2: maximize theta and phi

#include <iostream>
#include <fstream>
#include <cmath>
#include <joker-fdm/bvp1d.h>

using namespace std;

const double L = 50;
const double T = 30;
const double a = 0.92;
const double alpha = 10.0/3;
const double kappa_a = 0.01;
const double b = 18.7;
const double beta = 10;
const double theta_b1 = 0.3;
const double theta_b2 = 0.8;
const double c = 3e10;
// initial conditions
const double theta_init = 1;
const double phi_init = 1;

// minimal and maximal admissible controls
const double u_min = 0.01;
const double u_max = 0.5;

const int N = 500;  // the number of spatial grid points
const int M = 5000;  // the number of temporal grid points

const int iter_num = 5;  // number of iterations of the simple iterative method
const double u_init = u_min;  // initial guess

double u_left[M + 1], u_right[M + 1];  // values of the control
double J_values[iter_num];  // values of the objective functional

const char* output_file_name[2] = {"output_min.txt", "output_max.txt"};
const char* output_theta_phi_file_name[2] = {"output_theta_phi_min.txt",
    "output_theta_phi_max.txt"};
const char* output_J_file_name[2] = {"output_J_min.txt", "output_J_max.txt"};

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
    const double phi_b1 = pow(theta_b1, 4);
    const double phi_b2 = pow(theta_b2, 4);
    data.a[0][0] = a;
    data.a[1][0] = alpha;
    data.b[0][0] = data.b[0][1] = beta;
    data.w[0][0] = beta * theta_b1;
    data.w[0][1] = beta * theta_b2;
    // data.b[1][0] = ...;  data.b[1][1] = ...;
    // data.w[1][0] = ...;  data.w[1][1] = ...;
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

    for (int problem = 1; problem <= 2; ++problem) {
        // initialize the control
        for (int m = 0; m <= M; ++m)
            u_left[m] = u_right[m] = u_init;

        for (int iter = 1; iter <= iter_num; ++iter) {
            cout << "iteration " << iter << endl;
            // initial conditions
            for (int n = 0; n <= N; ++n) {
                sol[0](0, n) = theta_init;
                sol[1](0, n) = phi_init;
                for (int i = 0; i < 2; ++i)
                    sol_time[i](0, 0, n) = sol[i](0, n);
            }
            // solve the nonstationary problem
            for (int m = 1; m <= M; ++m) {
                cout << "m = " << m << " of " << M << endl;
                data.b[1][0] = u_left[m];
                data.b[1][1] = u_right[m];
                data.w[1][0] = u_left[m] * phi_b1;
                data.w[1][1] = u_right[m] * phi_b2;
                for (int n = 0; n <= N; ++n) {
                    data.g[0](0, n) = cvec[0] / tau * sol[0](0, n);
                    data.g[1](0, n) = cvec[1] / tau * sol[1](0, n);
                }
                SolveBVP1D(data, Parameters1D(), sol);
                for (int i = 0; i < 2; ++i) {
                    for (int n = 0; n <= N; ++n)
                        sol_time[i](m, 0, n) = sol[i](0, n);
                }
            }

            // update the control
            if (problem == 1) {
                for (int m = 0; m <= M; ++m) {
                    if (sol_time[1](m, 0, 0) - phi_b1 < 0)
                        u_left[m] = u_min;
                    else
                        u_left[m] = u_max;
                    if (sol_time[1](m, 0, N) - phi_b2 < 0)
                        u_right[m] = u_min;
                    else
                        u_right[m] = u_max;
                }
            }
            else {  // problem == 2
                for (int m = 0; m <= M; ++m) {
                    if (sol_time[1](m, 0, 0) - phi_b1 > 0)
                        u_left[m] = u_min;
                    else
                        u_left[m] = u_max;
                    if (sol_time[1](m, 0, N) - phi_b2 > 0)
                        u_right[m] = u_min;
                    else
                        u_right[m] = u_max;
                }
            }

            // calculate the objective functional
            double J = 0;
            for (int m = 0; m <= M; ++m) {
                double s = 0;
                for (int n = 0; n <= N; ++n) {
                    double term = pow(sol_time[0](m, 0, n), 2);
                    if (n == 0 || n == N)
                        term /= 2;
                    s += term;
                }
                if (m == 0 || m == M)
                    s /= 2;
                J += s;
            }
            J *= tau * grid.h[0];
            J_values[iter - 1] = J;
        }

        // output the optimal control
        ofstream fout(output_file_name[problem - 1]);
        for (int m = 0; m <= M; ++m)
            fout << m * tau << "  " << u_left[m] << "  " << u_right[m] << endl;
        // output theta and phi at the final time instant
        ofstream fout_theta_phi(output_theta_phi_file_name[problem - 1]);
        fout_theta_phi.precision(10);
        for (int n = 0; n <= N; ++n) {
            fout_theta_phi << n * grid.h[0] << "  " << sol_time[0](M, 0, n) << "  "
                << sol_time[1](M, 0, n) << endl;
        }
        // output values of the objective functional
        ofstream fout_J(output_J_file_name[problem - 1]);
        fout_J.precision(10);
        for (int iter = 0; iter < iter_num; ++iter)
            fout_J << J_values[iter] << endl;
    }  // for problem

    return 0;
}
