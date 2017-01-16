# A nonlinear problem in 1D
# -u''(x) + u^2(x) = sin(x) + sin^2(x) on (0, 2*pi)
# u(0) = 0, u(2*pi) = 0
# Exact solution: u(x) = sin(x)

clear all;
more off;
format long;

data.N = N = 1;
data.M = M = 1;
data.a = [1];
data.f = {@(x) x ^ 2};
data.df = {@(x) 2 * x};
data.b = [Inf, Inf];
data.w = [0, 0];
data.G = [];
L = [2 * pi];
Ltot = sum(L);

gfun = @(x) sin(x) + sin(x) ^ 2;
uexact = @(x) sin(x);

addpath("../../bvp1d");
K = [10];
for test = 0 : 3
  grid_info = get_grid_info(L, K);
  xgrid = linspace(0, L(1), K(1) + 1);
  data.g = arrayfun(gfun, xgrid);
  guess = zeros(1, grid_info.nodes);
  sol = solve_bvp1d(grid_info, data, guess, 1e-8, 20);
  sol_exact = arrayfun(uexact, xgrid);
  rms = sqrt(meansq(sol .- sol_exact));
  printf("h = %e, %e\n", grid_info.h);
  printf("rms = %e\n", rms);
  if (test > 0)
    printf("rms_old / rms = %f\n", rms_old / rms);
  endif
  rms_old = rms;
  K .*= 2;
endfor

figure
sol_exact = arrayfun(uexact, xgrid);
plot(xgrid, sol, xgrid, sol_exact);
xlabel("x");
ylabel("u");
