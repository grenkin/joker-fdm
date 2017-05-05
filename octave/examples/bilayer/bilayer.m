# The same problem as at https://github.com/grenkin/bilayer-fdm

clear all
more off
format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Physical data for two layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data for layer 1
layer_data(1) = struct(
  "L",     1,
  "omega", 0.2,
  "n",     1.8,
  "Nc",    0.001,
  "A",       0);
% data for the left boundary
boundary_data(1) = struct(
  "R",      0.3,
  "thetab",   1);

% data for layer 2
layer_data(2) = struct(
  "L",       1,
  "omega", 0.9,
  "n",       1,
  "Nc",   0.0001,
  "A",       0);
% data for the right boundary
boundary_data(2) = struct(
  "R",      0.2,
  "thetab", 0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param_func

for j = 1 : 2
  boundary_data(j).n = layer_data(j).n;
  layer(j) = calc_layer(layer_data(j));
  boundary(j) = calc_boundary(boundary_data(j));
endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data.N = N = 2;
data.M = M = 2;
data.a(1, :) = [layer(:).C];
data.a(2, :) = [layer(:).B];

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

data.G = [ Inf ; calc_G(layer_data(1).n, layer_data(2).n) ];

L = [layer(:).L];
K = 100 * [1, 1];

addpath("../../bvp1d");
grid_info = get_grid_info(L, K);
data.g = zeros(N, grid_info.nodes);
guess = [ boundary(1).thetab + ...
  (boundary(2).thetab - boundary(1).thetab) * (0 : (grid_info.nodes - 1)) / ...
  (grid_info.nodes - 1) ; zeros(1, grid_info.nodes) ];
tol = 1e-6;

# Calculate with taking into account refraction and reflection
sol = solve_bvp1d(grid_info, data, guess, tol);
theta = sol(1, :);
phi = sol(2, :);

# Calculate without taking into account refraction and reflection
data.G = [ Inf ; Inf ];
sol1 = solve_bvp1d(grid_info, data, guess, tol);
theta1 = sol1(1, :);
phi1 = sol1(2, :);

xgrid1 = (0 : K(1)) * grid_info.h(1);
xgrid2 = L(1) + (0 : K(2)) * grid_info.h(2);
xgrid = [xgrid1, xgrid2];
figure
plot(xgrid, theta, "r", "linewidth", 2, ...
  xgrid, theta1, "color", [0 .6 0], "linewidth", 2);
xlabel("x");
ylabel("theta");

figure
plot(xgrid1, phi(1 : (grid_info.K(1) + 1)), "r", "linewidth", 2, ...
  xgrid2, phi(grid_info.first_index(2) : grid_info.nodes), "r", "linewidth", 2, ...
  xgrid, phi1, "color", [0 .6 0], "linewidth", 2);
xlabel("x");
ylabel("phi");
