# A nonlinear conjugation problem in 1D
# -u''(x) + u^4 = sin(x) + sin^4(x) on (0, pi/2)
# -4u''(x) + u = 2 * sin(x/2 + pi/4) on (pi/2, 3*pi/2)
# -u'(0) + u(0) = -1, 4u'(3*pi/2) + u(3*pi/2) = -2
# u'(pi/2 - 0) = 4u'(pi/2 + 0)
# Exact solution: u(x) = sin(x), x in (0, pi/2);
#                 u(x) = sin(x/2 + pi/4), x in (pi/2, 3*pi/2)

clear all;
more off;
format long;

data.N = N = 1;
data.M = M = 2;
data.a = [1, 4];
data.f = {@(x) x ^ 4, @(x) x};
data.df = {@(x) 4 * x ^ 3, @(x) 1};
data.b = [1, 1];
data.w = [-1, -2];
data.G = Inf;
L = [pi / 2, pi];
Ltot = sum(L);

gfun = @(x) ifelse(x < pi / 2, sin(x) + sin(x) ^ 4, 2 * sin(x / 2 + pi / 4));
uexact = @(x) ifelse(x < pi / 2, sin(x), sin(x / 2 + pi / 4));

addpath("../../bvp1d");
K = [10, 10];
for test = 0 : 3
  grid_info = get_grid_info(L, K);
  xgrid = [linspace(0, L(1), K(1) + 1), linspace(L(1), Ltot, K(2) + 1)];
  data.g = arrayfun(gfun, xgrid);
  guess = zeros(1, grid_info.nodes);
  sol = solve_bvp1d(grid_info, data, guess, 1e-8, 10);
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
