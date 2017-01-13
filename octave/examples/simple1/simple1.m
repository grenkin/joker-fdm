# A simple linear problem in 1D
# -u''(x) = sin(x) on (0, pi)
# -u'(0) + u(0) = 0, u'(pi) + u(pi) = 0
# Exact solution: u(x) = 1 + sin(x)

clear all;
more off;
format long;

data.N = N = 1;
data.M = M = 1;
data.a = 1;
data.f = data.df = {@(x) 0};
data.b = [1, 1];
data.w = [0, 0];
data.G = [];
L = pi;

gfun = @sin;
uexact = @(x) 1 + sin(x);

addpath("../../bvp1d");
K = 10;
for test = 0 : 5
  grid_info = get_grid_info(L, K);
  data.g = arrayfun(gfun, linspace(0, L, grid_info.nodes));
  guess = zeros(1, grid_info.nodes);
  sol = solve_bvp1d(grid_info, data, guess, 1e-8, 1);
  sol_exact = arrayfun(uexact, linspace(0, L, grid_info.nodes));
  rms = meansq(sol .- sol_exact);
  printf("h = %e\n", grid_info.h(1));
  printf("rms = %e\n", rms);
  if (test > 0)
    printf("rms_old / rms = %f\n", rms_old / rms);
  endif
  rms_old = rms;
  K *= 2;
endfor

figure
xgrid = linspace(0, L, grid_info.nodes);
sol_exact = arrayfun(uexact, xgrid);
plot(xgrid, sol, xgrid, sol_exact);
xlabel("x");
ylabel("u");
