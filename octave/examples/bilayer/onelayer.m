# The same problem as at https://github.com/grenkin/bilayer-fdm

clear all;
more off;
format long;

% data for the layer
layer_data(1) = struct(
  "L",     1,
  "omega", 0.9,
  "n",     1,
  "Nc",    1/600,
  "A",     1);
% data for the left boundary
boundary_data(1) = struct(
  "R",      0.3,
  "thetab", 1);
% data for the right boundary
boundary_data(2) = struct(
  "R",      0.3,
  "thetab", 0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_func

layer(1) = calc_layer(layer_data(1));
for j = 1 : 2
  boundary_data(j).n = layer_data(1).n;
  boundary(j) = calc_boundary(boundary_data(j));
endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.N = N = 2;
data.M = M = 1;
data.a(1, 1) = layer(1).C;
data.a(2, 1) = layer(1).B;

fm =  {@(x) x^4,  @(x) -x;
       @(x) -x^4, @(x) x};
dfm = {@(x) 4*x^3,  @(x) -1;
       @(x) -4*x^3, @(x) 1};
data.f = data.df = cell(N, M, N);
for j = 1 : M
  for i = 1 : N
    for k = 1 : N
      data.f{i, j, k} = @(x) layer(j).K * fm{i, k}(x);
      data.df{i, j, k} = @(x) layer(j).K * dfm{i, k}(x);
    endfor
  endfor
endfor

data.b = [ Inf, Inf ; [boundary(:).tilde_gamma] ];
data.w = [ [boundary(:).thetab] ; [boundary(:).tilde_gamma] .* [boundary(:).thetab] .^ 4 ];

data.G = [ ];

L = layer(1).L;
K = 1000;

addpath("../../bvp1d");
grid_info = get_grid_info(L, K);
data.g = zeros(N, grid_info.nodes);
guess = zeros(N, grid_info.nodes);
tol = 1e-8;
sol = solve_bvp1d(grid_info, data, guess, tol);
theta = sol(1, :);
phi = sol(2, :);

figure
xgrid = [(0 : K(1)) * grid_info.h(1)];
plot(xgrid, theta, "r", "linewidth", 2);
xlabel("x");
ylabel("theta");
