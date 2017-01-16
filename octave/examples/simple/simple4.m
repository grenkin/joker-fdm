# A simple linear conjugation problem in 1D
# -u''(x) = -2 on (0, 1), -2u''(x) = -2 on (1, 3)
# u(0) = 0, 2u'(3) + u(3) = 11, u'(1 - 0) = 2u'(1 + 0), u(1 - 0) = u(1 + 0)
# Exact solution: u(x) = x^2 on (0, 1), u(x) = (x^2 + 1)/2 on (1, 3)

clear all;
more off;
format long;

data.N = N = 1;
data.M = M = 2;
data.a = [1, 2];
data.f = data.df = {@(x) 0, @(x) 0};
data.b = [Inf, 1];
data.w = [0, 11];
data.G = Inf;
L = [1, 2];
Ltot = sum(L);

gfun = @(x) -2;
uexact = @(x) ifelse(x < 1, x ^ 2, (x ^ 2 + 1) / 2);

addpath("../../bvp1d");
K = [10, 10];
for test = 0 : 5
  grid_info = get_grid_info(L, K);
  xgrid = [linspace(0, L(1), K(1) + 1), linspace(L(1), Ltot, K(2) + 1)];
  data.g = arrayfun(gfun, xgrid);
  guess = zeros(1, grid_info.nodes);
  sol = solve_bvp1d(grid_info, data, guess, 1e-8, 1);
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
