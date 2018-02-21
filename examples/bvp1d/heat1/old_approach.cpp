#include <stdio.h>
#include "umfpack.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <math.h>
#include <iomanip>

const double L = 0.5;  // length of the layer
const double T = 6e-8;  // length of the time interval

// parameters of air at 400 C and normal atmospheric pressure
const double k = 0.0521;  // thermal conductivity
const double rho = 0.524;  // density
const double c_p = 1068;  // heat capacity
const double n = 1.0;  // refractive index

// parameters of the boundary
const double htr = 10;  // heat transfer coefficient
const double eps = 0.8;  // emissivity of the boundary

const double Tmax = 773;  // maximum temperature in the unnormalized model
const double sigma = 5.67e-8;  // Stefan-Boltzmann constant
const double a = k / (rho * c_p);
const double b = 4 * sigma * pow(n, 2) * pow(Tmax, 3) / (rho * c_p);

const double kappa = 10;  // extinction coefficient
const double kappa_a = 1;  // absorption coefficient
const double alpha = 1 / (3 * kappa);

const double beta = htr / (rho * c_p);
const double gamma = eps / (2 * (2 - eps));

const double c = 3e8;  // the speed of light
const double mu = 1 / c;

const double theta_b = 0.4;

const int N = 1000;  // the number of spatial grid points
const int tnum = 10000;  // the number of temporal grid points

// coefficients in the initial conditions
double A = 0.3;
double B = 0.3;

//initial functions

double zeta_init (double x)
{
    return A * (1 - cos(2 * M_PI * x / L));
}

double xi_init (double x)
{
    return B * (1 - cos(2 * M_PI * x / L));
}

// parameters of the iterative method
const int numiter_stat = 50;
const int numiter = 20;
const double eps_newton = 1e-13;

const int neq = 2*(N+1), neq_phi = N+1;
const int MAX_EQ_NONZEROS = 4;
const int MAX_NONZEROS = neq*MAX_EQ_NONZEROS;
const int MAX_COL_NONZEROS = 4;

typedef enum { THETA, PHI } var_t;
double theta[tnum+1][N+1], phi[tnum+1][N+1], phi_tilde[tnum+1][N+1];

double thetanew[N+1], thetaold[N+1], phinew[N+1];

int Ap[neq+1];
int Ai[MAX_NONZEROS];
double Ax[MAX_NONZEROS];
double rhs[neq];
double x[neq];

double mat[neq][MAX_EQ_NONZEROS];
int mat_col[neq][MAX_EQ_NONZEROS];
int mat_cnt[neq];

double mat1[neq][MAX_COL_NONZEROS];
int mat1_row[neq][MAX_COL_NONZEROS];
int mat1_cnt[neq];

int get_index(var_t v, int i)
{
    if (i < 0 || i > N) {
        std::cout << "\nBAD INDEX!!!\n";
        exit(1);
    }
    return (v == THETA  ?  i  :  N + 1 + i);
}

void fill_mat(var_t v1, int i1, var_t v2, int i2, double num)
{
    int eq = get_index(v1, i1);
    int c = get_index(v2, i2);
    mat[eq][mat_cnt[eq]] = num;
    mat_col[eq][mat_cnt[eq]] = c;
    ++mat_cnt[eq];
}

void fill_mat_phi(int i1, int i2, double num)
{
    fill_mat(THETA, i1, THETA, i2, num);
}

void fill_rhs(var_t v, int i, double num)
{
    int eq = get_index(v, i);
    rhs[eq] = num;
}

void fill_rhs_phi(int i, double num)
{
    fill_rhs(THETA, i, num);
}

void init_mat()
{
    for (int i = 0; i < neq; ++i)
        mat_cnt[i] = 0;
}

void convert_to_mat1(int n)
{
    for (int i = 0; i < n; ++i)
        mat1_cnt[i] = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < mat_cnt[i]; ++j) {
            int col = mat_col[i][j];
            mat1[col][mat1_cnt[col]] = mat[i][j];
            mat1_row[col][mat1_cnt[col]] = i;
            ++mat1_cnt[col];
        }
    }
}

void convert_to_comp_col(int n)
{
    int k = 0;
    for (int col = 0; col < n; ++col) {
        Ap[col] = k;
        for (int j = 0; j < mat1_cnt[col]; ++j) {
            Ax[k] = mat1[col][j];
            Ai[k] = mat1_row[col][j];
            ++k;
        }
    }
    Ap[n] = k;
}

double zetainit(double x)
{
    return A*( -cos(2*M_PI*x/L) + 1 );
}

double xiinit(double x)
{
    return B*( -cos(2*M_PI*x/L) + 1 );
}

int main (void)
{
    double h = L/N;
    double tau = T/tnum;

    double *null = (double *) NULL ;
    void *Symbolic, *Numeric ;

    double phi_b = pow(theta_b, 4);

    for (int i = 0; i <= N; ++i) {
        theta[0][i] = theta_b + zetainit(i*h);
        phi[0][i] = phi_b + xiinit(i*h);
    }

    for (int m = 1; m <= tnum; ++m) {
        std::cout << "m = " << m << "     ";

        for (int i = 0; i <= N; ++i)
            thetaold[i] = theta[m-1][i];

        for (int iter = 1; iter <= numiter; ++iter) {
            std::cout << iter << " ";
            init_mat();
            for (int i = 0; i <= N; ++i) {
                if (i == 0 || i == N) {
                    double betaval, gammaval, thetabval, phibval;
                    betaval = beta;
                    gammaval = gamma;
                    thetabval = theta_b;
                    phibval = phi_b;

                    int i1 = (i == 0  ?  1  :  N - 1);
                    fill_mat(THETA, i, THETA, i,
                        a/h + betaval + b*kappa_a*h/2*4*pow(thetaold[i], 3) + h/(2*tau));
                    fill_mat(THETA, i, THETA, i1,
                        -a/h);
                    fill_mat(THETA, i, PHI, i,
                        -b*kappa_a*h/2);
                    fill_rhs(THETA, i,
                        betaval*thetabval + b*kappa_a*h/2*3*pow(thetaold[i], 4) + h/(2*tau)*theta[m-1][i]);

                    fill_mat(PHI, i, PHI, i,
                        alpha/h + gammaval + kappa_a*h/2 + mu*h/(2*tau));
                    fill_mat(PHI, i, PHI, i1,
                        -alpha/h);
                    fill_mat(PHI, i, THETA, i,
                        -kappa_a*h/2*4*pow(thetaold[i], 3));
                    fill_rhs(PHI, i,
                        gammaval*phibval - kappa_a*h/2*3*pow(thetaold[i], 4) + mu*h/(2*tau)*phi[m-1][i]);
                }
                else { // 0 < i < N
                    fill_mat(THETA, i, THETA, i-1,
                        -a/(h*h));
                    fill_mat(THETA, i, THETA, i,
                        1/tau + 2*a/(h*h) + b*kappa_a*4*pow(thetaold[i], 3));
                    fill_mat(THETA, i, THETA, i+1,
                        -a/(h*h));
                    fill_mat(THETA, i, PHI, i,
                        -b*kappa_a);
                    fill_rhs(THETA, i,
                        1/tau*theta[m-1][i] + b*kappa_a*3*pow(thetaold[i], 4));

                    fill_mat(PHI, i, PHI, i-1,
                        -alpha/(h*h));
                    fill_mat(PHI, i, PHI, i,
                        mu/tau + 2*alpha/(h*h) + kappa_a);
                    fill_mat(PHI, i, PHI, i+1,
                        -alpha/(h*h));
                    fill_mat(PHI, i, THETA, i,
                        -kappa_a*4*pow(thetaold[i], 3));
                    fill_rhs(PHI, i,
                        mu/tau*phi[m-1][i] - kappa_a*3*pow(thetaold[i], 4));
                }
            }

            convert_to_mat1(neq);
            convert_to_comp_col(neq);

            (void) umfpack_di_symbolic (neq, neq, Ap, Ai, Ax, &Symbolic, null, null) ;
            (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
            umfpack_di_free_symbolic (&Symbolic) ;
            (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, rhs, Numeric, null, null) ;
            umfpack_di_free_numeric (&Numeric) ;

            for (int i = 0; i <= N; ++i) {
                int k = get_index(THETA, i);
                thetanew[i] = x[k];
                k = get_index(PHI, i);
                phinew[i] = x[k];
            }
            double maxdiff = 0.0;
            for (int i = 0; i <= N; ++i)
                maxdiff = fmax(maxdiff, fabs(thetanew[i] - thetaold[i]));
            if (maxdiff < eps_newton)
                break;
            for (int i = 0; i <= N; ++i)
                thetaold[i] = thetanew[i];
        } //for iter
        std::cout << "\n";
        for (int i = 0; i <= N; ++i) {
            theta[m][i] = thetanew[i];
            phi[m][i] = phinew[i];
        }
    } // for m

    // Solve the stationary equation for phi
    for (int m = 0; m <= tnum; ++m) {
        std::cout << "m = " << m << "     " << std::endl;
        init_mat();
        for (int i = 0; i <= N; ++i) {
            if (i == 0 || i == N) {
                double gammaval, phibval;
                gammaval = gamma;
                phibval = phi_b;
                int i1 = (i == 0  ?  1  :  N - 1);
                fill_mat_phi(i, i,
                    alpha/h + gammaval + kappa_a*h/2);
                fill_mat_phi(i, i1,
                    -alpha/h);
                fill_rhs_phi(i,
                    gammaval*phibval + kappa_a*h/2*pow(theta[m][i], 4));
            }
            else { // 0 < i < N
                fill_mat_phi(i, i-1,
                    -alpha/(h*h));
                fill_mat_phi(i, i,
                    2*alpha/(h*h) + kappa_a);
                fill_mat_phi(i, i+1,
                    -alpha/(h*h));
                fill_rhs_phi(i,
                    kappa_a * pow(theta[m][i], 4));
            }
        }

        convert_to_mat1(neq_phi);
        convert_to_comp_col(neq_phi);

        (void) umfpack_di_symbolic (neq_phi, neq_phi, Ap, Ai, Ax, &Symbolic, null, null) ;
        (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
        umfpack_di_free_symbolic (&Symbolic) ;
        (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, rhs, Numeric, null, null) ;
        umfpack_di_free_numeric (&Numeric) ;

        for (int i = 0; i <= N; ++i) {
            int k = get_index(THETA, i);
            phi_tilde[m][i] = x[k];
        }
    } // for m

    // calculate the norm of the difference
    std::vector<double> norm_sq(tnum + 1), r(tnum + 1);
    for (int m = 0; m <= tnum; ++m) {
        // use the trapezoid formula
        double integral =
            0.5 * pow(phi[m][0] - phi_tilde[m][0], 2)
            + 0.5 * pow(phi[m][N] - phi_tilde[m][N], 2);
        for (int n = 1; n <= N - 1; ++n)
            integral += pow(phi[m][n] - phi_tilde[m][n], 2);
        integral *= h;
        norm_sq[m] = integral;
        r[m] = sqrt(integral / L);
    }

    // output the result
    std::ofstream fmid("output_mid.txt");
    fmid.precision(20);
    int imid = N/2;
    for (int m = 0; m <= tnum; ++m) {
        fmid << m*tau << "   " << theta[m][imid] << "   " << phi[m][imid]
            << "   " << phi_tilde[m][imid] << std::endl;
    }
    std::ofstream fdiff("output_diff.txt");
    fdiff.precision(20);
    for (int m = 0; m <= tnum; ++m) {
        fdiff << m * tau << "  " << norm_sq[m] << "  " << r[m] << std::endl;
    }

    std::cout << "\nOK";

    return (0) ;
}
